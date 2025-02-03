#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:15:07 2022

@author: elia
"""

#function that performs a model selection procedure with k-fold cross validation
#and returns the trained model
import os
import json
import random
import numpy as np

from ML_resources.adapt_modelname import adapt_modelname
from ML_resources.get_model import get_model
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV

#controlling randomness
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value) #fixed ambient variable value
random.seed(seed_value) #fixed random module seed
np.random.seed(seed_value) #fixed numpy's random module seed


def model_train(x_train_scaled, y_train, modelname, seed_value, resources_path, 
                randomsearch = False,
                k = 10, num_thr = 2, optimize_for = "f1"): 
    
    """
    Trains a machine learning model using cross-validation.

    Args:
        x_train_scaled (array-like): Scaled input features for training.
        y_train (array-like): Target labels for training.
        modelname (str): Name of the model to train.
        seed_value (int): Seed value for reproducibility.
        resources_path (str): Path to the resources directory.
        randomsearch (bool, optional): Flag to indicate whether to perform random search. Defaults to False.
        k (int, optional): Number of cross-validation folds. Defaults to 10.
        num_thr (int, optional): Number of threads to use. Defaults to 2.
        optimize_for (str, optional): Metric to optimize during model selection. Defaults to "f1".

    Returns:
        tuple: A tuple containing the trained classifier, means, and standard deviations of the test scores.

    """
    
    modelname = adapt_modelname(modelname)
    
    if randomsearch:
        param_filename = os.path.join(resources_path,"RS{}_cv_params.txt".format(modelname))
    else:
        param_filename = os.path.join(resources_path,"{}_cv_params.txt".format(modelname))
    
    #get scikit-learn model object 
    model = get_model(modelname)
    
    with open(param_filename) as f:
        model_parameters = json.load(f)
    
    #upping WT class importance with different class weights
    model_parameters["class_weight"] = [{-2:1, -1: 1, 0: 1.8, 1: 1.6, 2: 1}]

    if modelname.lower() in ["svc", "svm"]:
        model_parameters["probability"] = [True]
    
    if randomsearch:
        classifier = RandomizedSearchCV(estimator = model, param_grid = model_parameters,
                                        cv = k, n_iter = 100,
                                        scoring = optimize_for, n_jobs = num_thr, verbose = 1)
    else:
        classifier = GridSearchCV(estimator = model, param_grid = model_parameters, cv = k, 
                                  scoring = optimize_for, n_jobs = num_thr, verbose = 1, error_score="raise")
        
    
    #train classifier
    print("Selecting best {} model".format(modelname))
    classifier.fit(x_train_scaled, y_train)
    
    print("Best Parameters: ", classifier.best_params_)
    #gathering the training results
    means = classifier.cv_results_["mean_test_score"]
    stds = classifier.cv_results_["std_test_score"]
    
    #ordering the results from worse to better
    sorted_pairs = sorted(zip(means, stds))
    tuples = zip(*sorted_pairs)
    means, stds = [list(tuple) for tuple in tuples]
    
    return classifier, means, stds
