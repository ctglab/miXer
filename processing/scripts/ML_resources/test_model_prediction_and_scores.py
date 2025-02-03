#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 20:28:58 2022

@author: elia
"""

#function that calculates metrics for train set, nosegdup and segdup
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
import numpy as np
from ML_resources.substitute_values import substitute_values

#TODO lentissimo, parallelizzare? 
def test_on_dataset(model, X, Y_true, avg, is_nsd_sd=False):
    """
    Tests a given model on a dataset by performing the following steps:
    1. Uses the model to make predictions and prediction probabilities for the input data (X).
    2. Calculates accuracy, precision, recall, and F1 score using the predicted labels and true labels (Y_true).

    Parameters:
    - model (object): The trained model used for prediction.
    - X (array-like): The input data to make predictions on.
    - Y_true (array-like): The true labels corresponding to the input data.
    - avg (str): The averaging method used for calculating precision, recall, and F1 score.
    - is_nsd_sd (bool): Indicates whether the dataset is NSD (non-standard deviation) or SD (standard deviation). Default is False.

    Returns:
    - result (list): A list containing the following elements:
        - prediction (array-like): The predicted labels for the input data.
        - allClasses_prediction_probabilities (array-like): The prediction probabilities for all classes.
        - accuracy (float): The accuracy score of the model.
        - precision (float): The precision score of the model.
        - recall (float): The recall score of the model.
        - f1 (float): The F1 score of the model.
    """
    
    #use model to get predictions and predictions probabilities
    prediction = model.predict(X)
    allClasses_prediction_probabilities = model.predict_proba(X)
    
    #changing prediction values if dataset is NSD or SD
    if is_nsd_sd:
        prediction = substitute_values(prediction)
        # Calculate class -1 probability as the sum of columns 0 and 1 (-2 and -1 pred proba)
        a = np.sum(allClasses_prediction_probabilities[:, [0, 1]], axis=1)
        # Extract class 0 probability as column 2
        b = allClasses_prediction_probabilities[:, 2]
        # Calculate class 1 probability as the sum of columns 3 and 4 (1 and 2+ pred proba)
        c = np.sum(allClasses_prediction_probabilities[:, [3, 4]], axis=1)
        # Create the new numpy array with columns [a, b, c]
        allClasses_prediction_probabilities = np.column_stack((a, b, c))
        
    #calculate accuracy, precision, recall and f1 score
    accuracy = accuracy_score(Y_true, prediction)
    precision = precision_score(Y_true, prediction, average = avg)
    recall = recall_score(Y_true, prediction, average = avg)
    f1 = f1_score(Y_true, prediction, average = avg)
    
    return [prediction, allClasses_prediction_probabilities, accuracy, precision, recall, f1]

def test_model_prediction_and_scores(model, X_train_golden_scaled, X_test_golden_scaled,
                                     X_nosegdup_scaled, X_segdup_scaled, 
                                     Y_train_golden, Y_test_golden, Y_nosegdup, Y_segdup,
                                     avg = "binary"):
    """
    Tests a given model on multiple datasets and returns the predictions and scores for each dataset.
    
    Parameters:
    - model: The trained model to be tested.
    - X_train_golden_scaled: The scaled input data for the golden training dataset.
    - X_test_golden_scaled: The scaled input data for the golden test dataset.
    - X_nosegdup_scaled: The scaled input data for the no-segdup dataset.
    - X_segdup_scaled: The scaled input data for the segdup dataset.
    - Y_train_golden: The true labels for the golden training dataset.
    - Y_test_golden: The true labels for the golden test dataset.
    - Y_nosegdup: The true labels for the no-segdup dataset.
    - Y_segdup: The true labels for the segdup dataset.
    - avg (str): The averaging method used for calculating precision, recall, and F1 score.
    
    Returns:
    - results (tuple): A tuple containing the following elements:
        - gld_train (list): A list of results for the golden training dataset obtained from `test_on_dataset`.
        - gld_test (list): A list of results for the golden test dataset obtained from `test_on_dataset`.
        - nsd (list): A list of results for the no-segdup dataset obtained from `test_on_dataset`.
        - sd (list): A list of results for the segdup dataset obtained from `test_on_dataset`.
    """
    
    gld_train = test_on_dataset(model, X_train_golden_scaled, Y_train_golden, avg, is_nsd_sd=False)
    print("XLR train done")
    gld_test = test_on_dataset(model, X_test_golden_scaled, Y_test_golden, avg, is_nsd_sd=False)
    print("XLR test done")
    #print("Predicting NSD")
    nsd = test_on_dataset(model, X_nosegdup_scaled, Y_nosegdup, avg, is_nsd_sd=True)
    print("NSD done")
    #print("Predicting SD")
    sd = test_on_dataset(model, X_segdup_scaled, Y_segdup, avg, is_nsd_sd=True)
    print("SD done")
    
    return gld_train, gld_test, nsd, sd
    