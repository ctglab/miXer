#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 12:24:07 2022

@author: elia
"""

#main script for miXer model training
import os
import sys
import argparse
import random
import math
import threading
import numpy as np
import pandas as pd
from joblib import load
import warnings

from ML_resources.str_to_bool import str_to_bool
from ML_resources.load_datasets import load_datasets
from ML_resources.dataset_setup import dataset_setup
from ML_resources.get_scaler_name import get_scaler_name
from ML_resources.get_trained_model_and_scaler import get_trained_model_and_scaler
from ML_resources.test_models import test_model
from ML_resources.split import split
from ML_resources.parallel_predictions import parallel_predictions
from ML_resources.save_run_params import save_run_params

#absolute path of miXer home directory
training_main_dir = os.path.dirname(os.path.abspath(__file__))
#absolute path of project directory
project_dir = os.path.abspath(os.path.join(training_main_dir, os.pardir))
#resources directory
resources_dir = os.path.join(training_main_dir, "ML_resources")
#ML model output directory
default_main_output_dir = os.path.join(training_main_dir, "training_output")
#Use case test output directory
default_usecase_output_folder = os.path.join(training_main_dir, "usecase_output")
sys.path.append(resources_dir)

DEFAULT_NUMTHREADS = 3
DEFAULT_TRAIN_SAMPLES = 13000
DEFAULT_TEST_PERCENTAGE = 0.2
#default_sim_min2_data = None #TODO eliminare dopo inclusione con merged golden
DEFAULT_FORCE_MNORM = True
DEFAULT_SKIP_TESTED = False
DEFAULT_SAVE_CHRX = True

DEFAULT_MTT = ["SVC"]
DEFAULT_USE_NOISE = True
DEFAULT_MU = 0
DEFAULT_SIGMA = 0.05
DEFAULT_KFOLD = 10
DEFAULT_OPTIMIZE_FOR = "f1_macro"
DEFAULT_AVERAGING = "macro"
DEFAULT_SCALERTYPE = "RobustScaler"
DEFAULT_FORCETRAINING = False
DEFAULT_SKIP_CHRXTEST = False
DEFAULT_VERBOSE_LVL = 1

#controlling randomness
SEED_VALUE = 42
os.environ['PYTHONHASHSEED']=str(SEED_VALUE) #fixed ambient variable value
random.seed(SEED_VALUE) #fixed random module seed
np.random.seed(SEED_VALUE) #fixed numpy's random module seed

ap = argparse.ArgumentParser()
################ Required arguments
#TODO setup default dataset directory pointing to merged dataset
ap.add_argument("-ddir", "--dataset_directory", required = False, help = "Directory of training dataset folder.")
ap.add_argument("-mdir", "--model_directory",required = False, default = None, help = "Directory in which to find the trained model and scaler.\nBoth must be present.")
################ To be set for most runs
ap.add_argument("-tst", "--test_directory", nargs = "*", required = False, help = "Test data directory. Multiple directories allowed.")
ap.add_argument("-tst_exp_tags", "--tst_experiment_tags", nargs = "*", default = None, help = "Experiment(s) name. Must have the same number of experiment tags and -tst folders. Default = tst_folder_name + '_TestX' with X from 0 to #test folders.")
ap.add_argument("-mtag", "--model_tag", default = None, help = "Model tag. Default = Training folder name.") #use an already present model tag and the script will reload the pretrained model (to force retraing, use -ft True argument)
ap.add_argument("-ntr", "--num_threads", default = DEFAULT_NUMTHREADS, required = False, help = "Number of threads to use for cross validation of models. Default = {}".format(DEFAULT_NUMTHREADS))
ap.add_argument("-tsamp", "--train_samples", default = DEFAULT_TRAIN_SAMPLES, required = False, help = "Specify the number of exons to use for training. Overrides -tperc argument. Default = {}".format(DEFAULT_TRAIN_SAMPLES) )
ap.add_argument("-tperc", "--test_percentage", default = DEFAULT_TEST_PERCENTAGE, required = False, help = "Fraction of golden set to be used for testing. Default = {:.2f}%%".format(DEFAULT_TEST_PERCENTAGE*100))
ap.add_argument("-savex", "--save_prepared_chrX_data", default = DEFAULT_SAVE_CHRX, required = False, help = "Whether to save the ready to use ChrX datasets. Default = {}".format(DEFAULT_SAVE_CHRX))
ap.add_argument("-savex_folder", "--save_folder_chrX_data", default = None, required = False, help = "Prepared ChrX data save folder.")
#ap.add_argument("-min2", "--sim_min2_data", default = default_sim_min2_data, help = "Path to simulated double deletions dataset. Default = {} (won't be used to train the ML model).".format(default_sim_min2_data))
#TODO rimuovere force_mnorm perché tanto ora è di default
ap.add_argument("-force_mnorm", "--force_test_median_normalization", default = DEFAULT_FORCE_MNORM, help = "Whether to force the median normalization of test samples. Default = {}".format(DEFAULT_FORCE_MNORM))
ap.add_argument("-output", "--test_output_folder_path", default = default_usecase_output_folder, required = False, help = "Path to test data's output folder. Default = {}".format(default_usecase_output_folder))
ap.add_argument("-skipTested", "--skip_sample_if_tested", default = DEFAULT_SKIP_TESTED, required = False, help = "Whether to skip the prediction on a sample. Useful if redoing the prediction but not retraining the model. Default = {}".format(DEFAULT_SKIP_TESTED))
#TODO ora cerca in training_output, lasciamo così?
ap.add_argument("-mlout", "--ml_train_output_folder_path", default = default_main_output_dir, required = False, help = "Path to trained model root directory. Default = {}".format(default_usecase_output_folder))
ap.add_argument("-skTest", "--SKIP_CHRX_TEST", default = DEFAULT_SKIP_CHRXTEST, required = False, help = "Whether to skip the ChrX portion of model test. Useful when reusing an already trained model.\nDefault = {}".format(DEFAULT_SKIP_CHRXTEST))
ap.add_argument("-vrb", "--verbose_level", default = DEFAULT_VERBOSE_LVL, required = False, help = "Integer to select the verbosity of this script. 0 is max silence, TBD is fullly verbose.\nDefault = {}".format(DEFAULT_VERBOSE_LVL))
################ Change if you know what you are doing
ap.add_argument("-mtt", "--models_to_train", nargs = "*", default = DEFAULT_MTT, required = False, help = "List of trainable models. Default = {}.\nAvailable models: RF = Random Forest\tNN = Multilayer Perceptron\tSVC = Support Vector Classifier\nInput example: RF NN SVC".format(DEFAULT_MTT[0]))
ap.add_argument("-noise", "--noise", default = DEFAULT_USE_NOISE, required = False, help = "Whether to add noise to training data or not. Default: {}".format(DEFAULT_USE_NOISE))
ap.add_argument("-mu", "--noise_mean", default = DEFAULT_MU, required = False, help = "Mean value of added white noise. Default = {}".format(DEFAULT_MU))
ap.add_argument("-sigma", "--noise_variance", default = DEFAULT_SIGMA, required = False, help = "Variance of added white noise. Default = {}".format(DEFAULT_SIGMA))
ap.add_argument("-k", "--k_fold", default = DEFAULT_KFOLD, required = False, help = "K fold value for cross validation. Default = {}".format(DEFAULT_KFOLD))
ap.add_argument("-opt", "--optimize_for", default = DEFAULT_OPTIMIZE_FOR, required = False, help = "Optimization metric for cross validation. Default = {}".format(DEFAULT_OPTIMIZE_FOR))
ap.add_argument("-avg", "--averaging", default = DEFAULT_AVERAGING, required = False, help = "Averaging method for multiclass precision, recall and f1 score evaluation. Default = {}".format(DEFAULT_AVERAGING))
ap.add_argument("-sct", "--scaler_type", default = DEFAULT_SCALERTYPE, required = False, help = "Select scaler to use. Available scalers: [RobustScaler, StandardScaler, MinMaxScaler], Default = {}".format(DEFAULT_SCALERTYPE))
ap.add_argument("-ft", "--force_training", default = DEFAULT_FORCETRAINING, required = False, help = "Force the retraining of the models if already present. Default = {}".format(DEFAULT_FORCETRAINING))

args = vars(ap.parse_args())
dataset_dir = args["dataset_directory"]
if dataset_dir is None:
    raise ValueError("Must specify a chrX data folder for training or a trained model folder. Interrupting.")
num_thr = int(args["num_threads"])
split_fraction = float(args["test_percentage"])
train_samples = args["train_samples"]
mu = float(args["noise_mean"])
sigma = float(args["noise_variance"])
if train_samples is not None:
    train_samples = int(train_samples)
k = int(args["k_fold"])
models_to_train = args["models_to_train"]
optimize_for = args["optimize_for"]
averaging = args["averaging"]
na_paths = args["test_directory"]
usecase_output_folder = args["test_output_folder_path"]
main_output_dir = args["ml_train_output_folder_path"]
force_training = str_to_bool(args["force_training"], "Force training")
scaler_type = args["scaler_type"]
add_noise = str_to_bool(args["noise"], "Use noise")
expnames = args["tst_experiment_tags"]
#path_to_min2 = args["sim_min2_data"]
model_tag = args["model_tag"]
force_test_median_normalization = str_to_bool(args["force_test_median_normalization"], "Force Test samples median normalization")
SKIP_CHRX_TEST = str_to_bool(args["SKIP_CHRX_TEST"], "Skip model testing on chrX")
skip_tested = str_to_bool(args["skip_sample_if_tested"], "Skip prediction if file already present")
verbose = int(args["verbose_level"])
save_chrx_data = str_to_bool(args["save_prepared_chrX_data"], "Save ready to use ChrX data")
out_chrx_data = args["save_folder_chrX_data"]

#deal with no test samples specified (train model only case)
if na_paths is None:
    na_paths = []

if expnames is None:
    expnames = []
    for i in range(len(na_paths)):
        # Rimuovi eventuali separatori finali nel percorso
        curr_path =na_paths[i]
        while curr_path.endswith(os.path.sep):
            curr_path = curr_path[:-1]
        expnames.append(os.path.basename(curr_path) + "_Test{}".format(i) )

#Defining test output folder
if not os.path.isdir(usecase_output_folder):
    os.mkdir(usecase_output_folder)
    
columns = ["GC_content", "Length",  "NRC_poolNorm", "Class"]
training_columns = ["GC_content", "Length","NRC_poolNorm"]
disp_columns = ["GC", "Len", "NRC"]

#training
########## Training and chrX test phase
if dataset_dir is None:
    foldername = "ChrX_dataFolder_unspecified"
else:
    foldername = os.path.dirname(dataset_dir).split("/")[-1]
    
#model output folder
if model_tag is None:
    current_output_folder = os.path.join(main_output_dir, foldername)
    test_final_output_folders = []
    for item in expnames:
        test_final_output_folders.append(os.path.join(usecase_output_folder, item + "_" + foldername))

else:
    current_output_folder = os.path.join(main_output_dir, model_tag)
    test_final_output_folders = []
    for item in expnames:
        test_final_output_folders.append(os.path.join(usecase_output_folder, item + "_" + model_tag))

#check if current_output_folder exists, if it doesn't, create it
if not os.path.exists(current_output_folder):
    os.makedirs(current_output_folder)
#load the datasets
gld, nsd, sd = load_datasets(dataset_dir)
#TODO raise error if gold is not found (?)
#TODO raise warning if nsd and/or sd are not found
if gld is not None:
    gld = gld[columns].copy()
else:
    SKIP_CHRX_TEST = True
if nsd is not None:
    nsd = nsd[columns].copy()
else:
    SKIP_CHRX_TEST = True
if sd is not None:
    sd = sd[columns].copy()
else:
    SKIP_CHRX_TEST = True #can't do the full set of tests # TODO change in test function
    
#get class labels
labels = list(gld["Class"].unique())
labels.sort()
#preparing datasets: separating samples from labels, adding noise, creating noisy training data
if out_chrx_data is None:
    out_chrx_data = os.path.join(current_output_folder, "Prepared_ChrX_datasets")
    if not os.path.isdir(out_chrx_data):
        os.mkdir(out_chrx_data)
    msg = "Prepared ChrX datasets save dir not specified. Output will be saved in :{}".format(out_chrx_data)
    warnings.warn(msg, UserWarning)
#saving both the datasets with model predictions and clean
out_chrx_data_no_pred = os.path.join(out_chrx_data, "No_predictions")
if not os.path.isdir(out_chrx_data_no_pred):
    os.mkdir(out_chrx_data_no_pred)
prepared_datasets_dict = dataset_setup(gld, nsd, sd,
                                       train_samples, split_fraction,
                                       mu, sigma, add_noise, SEED_VALUE,
                                       skip_training = SKIP_CHRX_TEST,
                                       save_results = save_chrx_data,
                                       out_dir = out_chrx_data_no_pred
                                       )
save_run_params(dataset_dir,num_thr,
             split_fraction, train_samples, mu, sigma, k, models_to_train,
             optimize_for, averaging, na_paths, usecase_output_folder, current_output_folder,
             force_training, scaler_type, add_noise, expnames, model_tag,
             force_test_median_normalization, SKIP_CHRX_TEST, skip_tested, verbose,
             save_chrx_data, out_chrx_data, test_final_output_folders)

if verbose > 0:
    print("Training and testing models with datasets from {}".format(dataset_dir))
#train the selected models
if verbose > 0:
    print("No pretrained model directory specified.")
    if not force_training:
        print("Will check for presence of trained model.")
model_scaler_tuples = []
for modelname in models_to_train:
    modelname = modelname.upper()
    model_id = modelname.upper() #TODO serve per i file dei parametri della gridsearch, modificare? come?
    #defining output folder name
    if train_samples is None:
        model_foldername = modelname + "_TestFr_" + str(split_fraction) 
    else:
        model_foldername = modelname + "_TrainSamp_" + str(train_samples)
    model_foldername = model_foldername + "_Noise_" + str(add_noise)
    model_output_path = os.path.join(current_output_folder, model_foldername)
    modelname = modelname + get_scaler_name(scaler_type)
    scalerpath = os.path.join(model_output_path, model_id + get_scaler_name(scaler_type) + ".joblib")
    #if output folder does not exist, create it
    if not os.path.exists(model_output_path):
        os.makedirs(model_output_path)
    MODEL_SAVENAME = "{}_cv_{}.joblib".format(modelname, optimize_for)
    modelpath = os.path.join(model_output_path, MODEL_SAVENAME)
    #define output folder for training datasets with current model prediction
    out_chrx_data_with_pred = os.path.join(out_chrx_data, "PredictionsWith_" + model_foldername)
    if not os.path.isdir(out_chrx_data_with_pred):
        os.mkdir(out_chrx_data_with_pred)
        
    if add_noise:
        #trains the model with both clean data and noisy data together
        if verbose > 0:
            print("Training the model with both clean data and noisy data together")
        clf, scaler = get_trained_model_and_scaler(modelpath, scalerpath, scaler_type, resources_dir,
                                                   modelname, model_id, force_training, model_output_path,
                                                   num_thr, k, optimize_for,
                                                   prepared_datasets_dict["x_train_with_noise"], prepared_datasets_dict["y_train_with_noise"], 
                                                   SEED_VALUE)
        #testing the trained model
        if not SKIP_CHRX_TEST:
            #testing the trained model
            test_model(clf, model_output_path, disp_columns, training_columns,
                    model_id, modelname, labels,
                    scaler.transform(prepared_datasets_dict["x_train_with_noise"]), scaler.transform(prepared_datasets_dict["x_test_clean_only"]),
                    scaler.transform(prepared_datasets_dict["x_nsd_or"]), scaler.transform(prepared_datasets_dict["x_sd_or"]),
                    prepared_datasets_dict["y_train_with_noise"], prepared_datasets_dict["y_test_clean_only"], 
                    prepared_datasets_dict["y_nsd"], prepared_datasets_dict["y_sd"], 
                    averaging, out_chrx_data_with_pred)
        else:
            if verbose > 0:
                print("Skipping model testing on ChrX.")
    else:
        #train the model with only the clean data
        if verbose > 0:
            print("Training the model on clean data only")
        clf, scaler = get_trained_model_and_scaler(modelpath, scalerpath, scaler_type, resources_dir,
                                                   modelname, model_id, force_training, model_output_path,
                                                   num_thr, k, optimize_for,
                                                   prepared_datasets_dict["x_train_clean_only"], prepared_datasets_dict["y_train_clean_only"], 
                                                   SEED_VALUE)
        #testing the trained model
        if not SKIP_CHRX_TEST:
            test_model(clf, model_output_path, disp_columns, training_columns,
                    model_id, modelname, labels,
                    scaler.transform(prepared_datasets_dict["x_train_clean_only"]), scaler.transform(prepared_datasets_dict["x_test_clean_only"]),
                    scaler.transform(prepared_datasets_dict["x_nsd_or"]), scaler.transform(prepared_datasets_dict["x_sd_or"]),
                    prepared_datasets_dict["y_train_clean_only"], prepared_datasets_dict["y_test_clean_only"], 
                    prepared_datasets_dict["y_nsd"], prepared_datasets_dict["y_sd"], 
                    averaging, out_chrx_data_with_pred)
        else:
            if verbose > 0:
                print("Skipping model testing on ChrX.")
    model_scaler_tuples.append((model_id, modelname, clf,scaler))
    #saving file with run parameters
    
for item in model_scaler_tuples:
    curr_mid = item[0]
    curr_mname = item[1]
    curr_clf = item[2]
    curr_scaler = item[3]
    if verbose > 0:
        print("Making predictions on usecase samples.\nApplying median normalization of test samples NRC_poolNorm.")
        print("Model: {}".format(curr_mid))
    ################ Inference phase
    for i in range(len(na_paths)): #TODO cambiare nome a na_paths (target sample path?)
        na_path = na_paths[i]
        curr_out_folder = test_final_output_folders[i]
        #load ID_TARGET.tzt.gz dataset if test directory is specified
        t = []
        if na_path is not None:
            for parent, dirs, files in os.walk(na_path):
                t.extend(files)
                break
        #keep only target filenames
        t_names = [x for x in t if "TARGET.txt.gz" in x]
        if len(t_names) != 0:
            if verbose > 0:
                print("Samples from: {}".format(na_path))
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
                threads.append(threading.Thread(target=parallel_predictions, args= (na_path, item, add_noise,
                                                                                    curr_out_folder, foldername, curr_mname,
                                                                                    training_columns, curr_scaler, curr_clf, 
                                                                                    force_test_median_normalization,
                                                                                    train_samples, split_fraction, skip_tested),
                                                name = "WES samples single-exon {} prediction thread {}".format(curr_mname, THR) ) )
                THR = THR+1
            for t in threads:
                t.start()
            for t in threads:
                t.join()
