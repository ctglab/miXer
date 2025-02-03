#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:06:26 2022

@author: elia
"""


#functions that returns the trained (or retrained, or reloaded) ML model of choice
#then saves it in the correct folder

import os
import warnings
from joblib import load, dump

from ML_resources.get_scaler import get_scaler
from ML_resources.model_train import model_train
from ML_resources.plot_training_scores import plot_training_scores

def do_work(x_train, y_train, model_id, seed_value, resources_dir, scaler_path_to_file,
            scaler_type, model_path_to_file, modelname, model_output_folder_path,
            verbose=1, force_training=False, randomsearch=None, k=None, num_thr=None, optimize_for=None, scaler_only = False):
    
    """
    Performs the necessary steps for training a model and saving the results.

    Args:
        x_train (array-like): Training data features.
        y_train (array-like): Training data labels.
        model_id (str): Identifier for the model.
        seed_value (int): Seed value for random number generation.
        resources_dir (str): Path to the directory containing necessary resources.
        scaler_path_to_file (str): Path to save the scaler object.
        scaler_type (str): Type of scaler to use.
        model_path_to_file (str): Path to save the trained model.
        modelname (str): Name of the model.
        model_output_folder_path (str): Path to the output folder for saving model-related files.
        verbose (int, optional): Verbosity level. Defaults to 1.
        force_training (bool, optional): Flag to force retraining of the model. Defaults to False.
        randomsearch (dict, optional): Parameters for random search. Defaults to None.
        k (int, optional): Number of folds for cross-validation. Defaults to None.
        num_thr (int, optional): Number of threads for parallel processing. Defaults to None.
        optimize_for (str, optional): Metric to optimize for. Defaults to None.
        scaler_only (bool, optional): Whether to train only the scaler.

    Returns:
        tuple: Trained model and scaler objects if not scaler_only, otherwise scaler object only.
    """
    
    scaler = get_scaler(scaler_type)
    scaler.fit(x_train)
    dump(scaler, scaler_path_to_file)
    x_train_scaled = scaler.transform(x_train)
    
    if not scaler_only:
        if verbose > 1:
            print("Training the {} model.".format(modelname))
        
        clf, cv_means, cv_stds = model_train(x_train_scaled, y_train, model_id, seed_value, resources_dir,
                                             randomsearch=randomsearch, k=k, num_thr=num_thr,
                                             optimize_for=optimize_for)
        # Save the model
        dump(clf, os.path.join(model_path_to_file))
        # TODO: Save the best parameters
        if verbose > 1:
            print("Selected {} model parameters:".format(modelname))
            print(clf.best_params_)
        # Plot classifier performance during training
        plot_training_scores(cv_means, cv_stds, modelname, optimize_for, model_output_folder_path, color="blue", lw=2)
        
        return clf, scaler
    else:
        return scaler

def get_trained_model_and_scaler(model_path_to_file, scaler_path_to_file, scaler_type, resources_dir,
                                 modelname, model_id, force_training, model_output_folder_path,
                                 num_thr, k, optimize_for,
                                 x_train_or, y_train,
                                 seed_value,
                                 randomsearch=None,
                                 verbose=1):
    
    """
    Retrieves a trained model and scaler objects. If not present, performs training and saves the results.

    Args:
        model_path_to_file (str): Path to the trained model file.
        scaler_path_to_file (str): Path to the scaler file.
        scaler_type (str): Type of scaler to use.
        resources_dir (str): Path to the directory containing necessary resources.
        modelname (str): Name of the model.
        model_id (str): Identifier for the model.
        force_training (bool): Flag to force retraining of the model.
        model_output_folder_path (str): Path to the output folder for saving model-related files.
        num_thr (int): Number of threads for parallel processing.
        k (int): Number of folds for cross-validation.
        optimize_for (str): Metric to optimize for.
        x_train_or (array-like): Original training data features.
        y_train (array-like): Training data labels.
        randomsearch (dict): Parameters for random search.
        seed_value (int): Seed value for random number generation.
        verbose (int, optional): Verbosity level. Defaults to 1.

    Returns:
        tuple: Trained model and scaler objects.
    """
    
    # Check for already trained model
    if os.path.isfile(model_path_to_file) and not force_training:
        # Load model (for testing phase)
        clf = load(model_path_to_file)
        if verbose > 0:
            print("Trained {} model loaded".format(modelname))
        
        # Check if scaler is present
        if os.path.isfile(scaler_path_to_file):
            scaler = load(scaler_path_to_file)
            if verbose > 0:
                print("Trained scaler loaded.")
        else:
            # Scaler is not present, retrain it
            warnings.warn("Trained scaler not present. Trained {} model is present.\n"
                          "This should not happen, please check for {}.joblib file.\n"
                          "Training a new one now using the training data\n.".format(modelname, scaler_type))
            _, scaler = do_work(x_train_or, y_train, model_id, seed_value, resources_dir,
                                scaler_path_to_file, scaler_type, model_path_to_file, modelname,
                                model_output_folder_path, verbose=verbose, force_training=True,
                                randomsearch=randomsearch, k=k, num_thr=num_thr, optimize_for=optimize_for)[1]

        
    else:
        if verbose > 0:
            print("Trained {} model not present.".format(modelname))
            print("Training the scaler")
        
        clf, scaler = do_work(x_train_or, y_train, model_id, seed_value, resources_dir,
                              scaler_path_to_file, scaler_type, model_path_to_file, modelname,
                              model_output_folder_path, verbose=verbose, force_training=force_training,
                              randomsearch=randomsearch, k=k, num_thr=num_thr, optimize_for=optimize_for)
    
    # Save model hyperparameters to a file
    txtname = "{}_model_perf_report.txt".format(modelname)
    with open(os.path.join(model_output_folder_path, txtname), "w+") as f:
        f.write("{} model selected hyperparameters:\n".format(modelname))
        f.write(str(clf.best_params_))
    
    return clf, scaler


