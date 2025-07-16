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
# ap.add_argument("-tst", "--test_directory", nargs = "*", required = True, help = "Test data directory. Multiple directories allowed.")
ap.add_argument("-mdir", "--model_directory",required = False, default = DEFAULT_MODEL_DIR, help = "Directory in which to find the trained model and scaler.\nBoth must be present.")

# ap.add_argument("-tst_exp_tags", "--tst_experiment_tags", nargs = "*", default = None, help = "Experiment(s) name. Must have the same number of experiment tags and -tst folders. Default = tst_folder_name + '_TestX' with X from 0 to #test folders.")
ap.add_argument("-output", "--test_output_folder_path", default = default_usecase_output_folder, required = False, help = "Path to test data's output folder. Default = {}".format(default_usecase_output_folder))
ap.add_argument("-ntr", "--num_threads", default = DEFAULT_NUMTHREADS, required = False, help = "Number of threads to use for cross validation of models. Default = {}".format(DEFAULT_NUMTHREADS))
ap.add_argument("-force_mnorm", "--force_test_median_normalization", default = DEFAULT_FORCE_MNORM, help = "Whether to force the median normalization of test samples. Default = {}".format(DEFAULT_FORCE_MNORM))
ap.add_argument("-skipTested", "--skip_sample_if_tested", default = DEFAULT_SKIP_TESTED, required = False, help = "Whether to skip the prediction on a sample. Useful if redoing the prediction but not retraining the model. Default = {}".format(DEFAULT_SKIP_TESTED))
ap.add_argument("-vrb", "--verbose_level", default = DEFAULT_VERBOSE_LVL, required = False, help = "Integer to select the verbosity of this script. 0 is max silence, TBD is fullly verbose.\nDefault = {}".format(DEFAULT_VERBOSE_LVL))
ap.add_argument("-chrx_dataFolder", "--chrX_dataFolder", default = DEFAULT_CHRX_DATAFOLDER, required = False, help = "Used to specify the model training data source in the prediction folder name. Default = {}".format(DEFAULT_CHRX_DATAFOLDER))

args = vars(ap.parse_args())
with open(args['json'], 'r') as j:
    config = json.load(j)
PREPARED_SVM_DIR = os.path.join(config['main_outdir_host'], "*_datasets_testing_*")
model_directory = args["model_directory"]
num_thr = int(config["threads"])
logging.info(f"PREPARED svm dir is: {PREPARED_SVM_DIR}")
na_paths = glob.glob(PREPARED_SVM_DIR) #TODO: simplify
logging.info(f"Detected {len(na_paths)} samples from the given folder: these are {na_paths}")
usecase_output_folder = args["test_output_folder_path"]
expnames = config['exp_id']
force_test_median_normalization = str_to_bool(args["force_test_median_normalization"], "Force Test samples median normalization")
skip_tested = str_to_bool(args["skip_sample_if_tested"], "Skip prediction if file already present")
verbose = int(args["verbose_level"])
foldername = str(args["chrX_dataFolder"])

#Defining test output folder
if not os.path.isdir(usecase_output_folder):
    os.makedirs(usecase_output_folder)

if expnames is None:
    expnames = []
    for i in range(len(na_paths)):
        curr_path = na_paths[i]
        while curr_path.endswith(os.path.sep):
            curr_path = curr_path[:-1]
        expnames.append(os.path.basename(curr_path) + "_Test{}".format(i) )
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
    #Defining output folders
    test_final_output_folders = []
    for item in expnames:
        test_final_output_folders.append(os.path.join(usecase_output_folder, item + "_" + curr_mid))
    ################ Inference phase
    for i in range(len(na_paths)): #TODO cambiare nome a na_paths (target sample path?)
        na_path = na_paths[i]
        logging.debug(f"na path is {na_path}")
        curr_out_folder = test_final_output_folders[i]
        #load ID_TARGET.tzt.gz dataset if test directory is specified
        t = []
        if na_path is not None:
            for parent, dirs, files in os.walk(na_path):
                t.extend(files)
                break
        logging.debug(f"This would work on these files: {t}")
        #keep only target filenames
        t_names = [x for x in t if "TARGET.txt.gz" in x]
        if len(t_names) != 0:
            if verbose > 0:
                logging.info("Samples from: {}".format(na_path))
            #### parallelizing test samples single-exon predictions
            if num_thr > len(t_names):
                target_thrds = len(t_names)
            else:
                target_thrds = num_thr
            thread_chunk_size = math.floor(len(t_names)/ target_thrds)
            target_lists = split(t_names, thread_chunk_size)
            threads = []
            THR = 1
            for item in target_lists:
                threads.append(threading.Thread(target=parallel_predictions, args= (na_path, item, inferred_noise,
                                                                                    curr_out_folder, foldername, curr_mname,
                                                                                    training_columns, curr_scaler, curr_clf, 
                                                                                    force_test_median_normalization,
                                                                                    inferred_trainsamples, inferred_splitFraction, skip_tested),
                                                name = "WES samples single-exon {} prediction thread {}".format(curr_mname, THR) ) )
                THR = THR+1
            for t in threads:
                t.start()
            for t in threads:
                t.join()

















