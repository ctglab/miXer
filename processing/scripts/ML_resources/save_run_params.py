# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 21:22:43 2023

@author: ecero
"""

import datetime
import os

def save_run_params(dataset_dir,num_thr,
             split_fraction, train_samples, mu, sigma, k, models_to_train,
             optimize_for, averaging, na_paths, usecase_output_folder, main_output_dir,
             force_training, scaler_type, add_noise, expnames, model_tag,
             force_test_median_normalization, SKIP_CHRX_TEST, skip_tested, verbose,
             save_chrx_data, out_chrx_data, test_final_output_folders):
    """
    Save the run parameters and configurations to a text file.

    Parameters:
    - dataset_dir (str): Directory path to the dataset.
    - num_thr (int): Maximum number of threads specified.
    - split_fraction (float): Training dataset split fraction (Test set percentage over total).
    - train_samples (int): Number of training samples to use for the training set (overrides split fraction if specified).
    - mu (float): Mean of added noise.
    - sigma (float): Variance of added noise.
    - k (int): Number of folds for cross-validation.
    - models_to_train (str): ML Model(s) specified for training.
    - optimize_for (str): ML models optimization metric.
    - averaging (str): Performance calculation mode.
    - na_paths (str): Paths to test samples.
    - usecase_output_folder (str): Usecase output folder.
    - main_output_dir (str): Training results output directory.
    - force_training (bool): Force model training.
    - scaler_type (str): Scaler type.
    - add_noise (bool): Add noise flag.
    - expnames (str): Test sample(s) experiment name(s).
    - model_tag (str): Model tag.
    - force_test_median_normalization (bool): Force test sample median normalization of NRC_poolNorm.
    - SKIP_CHRX_TEST (bool): Skip ChrX model test.
    - skip_tested (bool): Skip already tested usecase samples.
    - verbose (int): Verbose level.
    - save_chrx_data (bool): Save prepared ChrX datasets.
    - out_chrx_data (str): Prepared ChrX dataset folder.
    - test_final_output_folders (str): List containing all the test output folders
    
    Returns:
    - None
    """
    #get current date
    # Get the current date and time
    current_datetime = datetime.datetime.now()
    # Format the date and time as "DD-MM-YYYY hour-minute-second"
    date = current_datetime.strftime("%d-%m-%Y %H-%M-%S")
    #open txt file
    f = open(os.path.join(main_output_dir, "Run_params_log_{}.txt".format(date)), "w")
    f.write("Dataset directory: {}\n".format(dataset_dir))
    f.write("Max number of threads specified: {}\n".format(str(num_thr)))
    f.write("Training dataset split fraction (Test set percentage over total): {}\n".format(str(split_fraction)))
    f.write("Number of training samples to use for training set (overrides split fraction if specified): {}\n".format(str(train_samples)))
    f.write("Add noise bool: {}\n".format(str(add_noise)))
    f.write("Mean of added noise (mu): {}\n".format(str(mu)))
    f.write("Variance of added noise (sigma): {}\n".format(str(sigma)))
    f.write("ML Model(s) specified for training: {}\n".format(" ".join(models_to_train)))
    f.write("ML models optimization metric: {}\n".format(optimize_for))
    f.write("Performance calculation mode: {}\n".format(averaging))
    f.write("Number of folds for cross-validation (k): {}\n".format(str(k)))
    f.write("Force model training: {}\n".format(str(force_training)))
    f.write("Scaler type: {}\n".format(scaler_type))
    f.write("Model tag: {}\n".format(model_tag))
    f.write("Skip ChrX model test: {}\n".format(str(SKIP_CHRX_TEST)))
    f.write("Save prepared ChrX datasets: {}\n".format(str(save_chrx_data)))
    f.write("Prepared ChrX dataset folder: {}\n".format(out_chrx_data))
    f.write("Training results output directory: {}\n".format(main_output_dir))
    f.write("Paths to test samples: {}\n".format(na_paths))
    f.write("Usecase output folder: {}\n".format(usecase_output_folder))
    f.write("Processed usecases output folders: {}\n".format(" ".join(test_final_output_folders)))
    f.write("Test sample(s) experiment name(s): {}\n".format(" ".join(expnames)))
    f.write("Force test sample median normalization of NRC_poolNorm: {}\n".format(str(force_test_median_normalization)))
    f.write("Skip already tested usecase samples: {}\n".format(str(skip_tested)))
    f.write("Verbose level: {}\n".format(verbose))
    #TODO when re-running, the script does not accept this do_it_again
    do_it_again = """Run it again? | python ML_training.py -ddir {}
                    -tst {} -tst_exp_tags {} -mtag {} -ntr {} -tsamp {} -tperc {}
                    -savex {} -savex_folder {} -force_mnorm {} -output {} -skipTested {}
                    -mlout {} -skTest {} -vrb {} -mtt {} -noise {} -mu {} -sigma {}
                    -k {} -opt {} -avg {} -sct {} -ft {}""".format(dataset_dir, " ".join(na_paths), " ".join(expnames),
                                model_tag, str(num_thr), str(train_samples), str(split_fraction),
                                str(save_chrx_data), out_chrx_data, str(force_test_median_normalization),
                                usecase_output_folder, str(skip_tested), main_output_dir, str(SKIP_CHRX_TEST),
                                str(verbose), " ".join(models_to_train), str(add_noise), str(mu), str(sigma),
                                str(k), optimize_for, averaging, scaler_type, str(force_training))
    f.write(do_it_again)
    f.close()
    return