#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 09:23:16 2023

@author: elia
"""

import os
import sys
import argparse
import random
import math
import threading
import glob
import numpy as np
import pandas as pd
import json
from joblib import load
import logging

# setup logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler('mixerDataset.log'),
        logging.StreamHandler(sys.stdout)
    ])

from ML_resources.str_to_bool import str_to_bool
#from ML_resources.get_scaler_name import get_scaler_name
#from ML_resources.get_trained_model_and_scaler import get_trained_model_and_scaler
from ML_resources.split import split
from ML_resources.parallel_predictions import parallel_predictions

#absolute path of miXer home directory
training_main_dir = os.path.dirname(os.path.abspath(__file__))
#absolute path of project directory
project_dir = os.path.abspath(os.path.join(training_main_dir, os.pardir))
#resources directory
resources_dir = os.path.join(training_main_dir, "ML_resources")
#Use case test output directory
default_usecase_output_folder = os.path.join(training_main_dir, "usecase_output")
sys.path.append(resources_dir)

#Default values
DEFAULT_NUMTHREADS = 3
DEFAULT_FORCE_MNORM = True
DEFAULT_SKIP_TESTED = True
DEFAULT_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "release_mdl", "SVC_TrainSamp_13000_Noise_True")
DEFAULT_VERBOSE_LVL = 1
DEFAULT_CHRX_DATAFOLDER = ""
training_columns = ["GC_content", "Length","NRC_poolNorm"]

#controlling randomness
SEED_VALUE = 42
os.environ['PYTHONHASHSEED']=str(SEED_VALUE) #fixed ambient variable value
random.seed(SEED_VALUE) #fixed random module seed
np.random.seed(SEED_VALUE) #fixed numpy's random module seed

ap = argparse.ArgumentParser()
ap.add_argument('-j', '--json', help="Path to the miXer json file", required=True)
ap.add_argument("-mdir", "--model_directory",required = False, default = DEFAULT_MODEL_DIR, help = "Directory in which to find the trained model and scaler.\nBoth must be present.")
ap.add_argument("-force_mnorm", "--force_test_median_normalization", default = DEFAULT_FORCE_MNORM, help = "Whether to force the median normalization of test samples. Default = {}".format(DEFAULT_FORCE_MNORM))
ap.add_argument("-skipTested", "--skip_sample_if_tested", default = DEFAULT_SKIP_TESTED, required = False, help = "Whether to skip the prediction on a sample. Useful if redoing the prediction but not retraining the model. Default = {}".format(DEFAULT_SKIP_TESTED))
ap.add_argument("-vrb", "--verbose_level", default = DEFAULT_VERBOSE_LVL, required = False, help = "Integer to select the verbosity of this script. 0 is max silence, TBD is fullly verbose.\nDefault = {}".format(DEFAULT_VERBOSE_LVL))
ap.add_argument("-chrx_dataFolder", "--chrX_dataFolder", default = DEFAULT_CHRX_DATAFOLDER, required = False, help = "Used to specify the model training data source in the prediction folder name. Default = {}".format(DEFAULT_CHRX_DATAFOLDER))

args = vars(ap.parse_args())
with open(args['json'], 'r') as j:
    config = json.load(j)
PREPARED_SVM_DIR = os.path.join(
    os.path.abspath(config['main_outdir_host']),
    config['exp_id'],
    "datasets_miXer")
model_directory = args["model_directory"]
num_thr = int(config["threads"])
logging.info(f"PREPARED svm dir is: {PREPARED_SVM_DIR}")
sample_files = sorted([f for f in os.listdir(PREPARED_SVM_DIR) if f.endswith("_miXer_data.tsv")]) if os.path.isdir(PREPARED_SVM_DIR) else []
logging.info(f"Detected {len(sample_files)} samples from the given folder: these are {sample_files}")
expnames = [f.replace("_miXer_data.tsv", "") for f in sample_files]
usecase_output_folder = os.path.join(
    os.path.abspath(config['main_outdir_host']),
    config['exp_id'])
force_test_median_normalization = str_to_bool(args["force_test_median_normalization"], "Force Test samples median normalization")
skip_tested = str_to_bool(args["skip_sample_if_tested"], "Skip prediction if file already present")
verbose = int(args["verbose_level"])
foldername = str(args["chrX_dataFolder"])

#Defining test output folder
if not os.path.isdir(usecase_output_folder):
    os.makedirs(usecase_output_folder)

logging.info(f"Generated these expnames: {expnames}")
#reloading the model in the specified folder
logging.info("Reloading trained model from folder {}".format(model_directory))
model_scaler_tuples = []
#Inferring model tag from model directory
splits = [x for x in os.path.normpath(model_directory).split(os.path.sep) if x != '']
inferred_model_folder = splits[-1]
inferred_model_tag = inferred_model_folder.split("_")[0]
inferred_noise = inferred_model_folder.split("_")[4]
inferred_trainsamples_or_splitFraction = float(inferred_model_folder.split("_")[2])
if inferred_trainsamples_or_splitFraction > 1:
    inferred_trainsamples = inferred_trainsamples_or_splitFraction
    inferred_splitFraction  =None
else:
    inferred_trainsamples = None
    inferred_splitFraction = inferred_trainsamples_or_splitFraction

joblib_files = [x for x in os.listdir(model_directory) if ".joblib" in x]
scalerfile = [x for x in joblib_files if "cv" not in x][0]
modelfile = [x for x in joblib_files if "cv" in x][0]
inferred_modelname = modelfile.split("_")[0]
logging.info("Found {} model in folder {}\nModel file: {} | Scaler file: {}".format(inferred_model_tag,
                                                                             inferred_model_folder,
                                                                             modelfile, scalerfile))
clf = load(os.path.join(model_directory, modelfile))
scaler = load(os.path.join(model_directory, scalerfile))
model_scaler_tuples.append((inferred_model_tag, inferred_modelname, clf, scaler))


#Inference step with model-scaler reloaded tuples
for item in model_scaler_tuples:
    curr_mid, curr_mname, curr_clf, curr_scaler = item
    if verbose > 0:
        logging.debug("Making predictions on usecase samples.\nApplying median normalization of test samples NRC_poolNorm.")
        logging.debug("Model: {}".format(curr_mid))
    ################ Inference phase
    if len(sample_files) != 0:
        logging.info("Samples from: {}".format(PREPARED_SVM_DIR))
        #### parallelizing test samples predictions
        if num_thr > len(sample_files):
            target_thrds = len(sample_files)
        else:
            target_thrds = num_thr
        thread_chunk_size = math.floor(len(sample_files) / target_thrds)
        target_lists = split(sample_files, thread_chunk_size)
        threads = []
        THR = 1
        for chunk in target_lists:
            threads.append(threading.Thread(target=parallel_predictions, args=(PREPARED_SVM_DIR, chunk, inferred_noise,
                                                                                usecase_output_folder, curr_mid, foldername, curr_mname,
                                                                                training_columns, curr_scaler, curr_clf,
                                                                                force_test_median_normalization,
                                                                                inferred_trainsamples, inferred_splitFraction, skip_tested),
                                            name="WES samples {} prediction thread {}".format(curr_mname, THR)))
            THR += 1
        for t in threads:
            t.start()
        for t in threads:
            t.join()

















