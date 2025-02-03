#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 23:04:19 2022

@author: elia
"""

#function that calculates all the scores for a single prediction-true tuple 
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

def training_model_scores(y_pred, y_true, avg="binary"):
    """
    Calculate various evaluation scores for a classification model.

    Args:
        y_pred (array-like): Predicted labels.
        y_true (array-like): True labels.
        avg (str, optional): The averaging strategy to use for precision, recall, and F1 score.
            Supported options are 'binary' (default), 'micro', 'macro', and 'weighted'.

    Returns:
        list: A list containing the following evaluation scores:
            - model_accuracy: Accuracy score.
            - model_precision: Precision score.
            - model_recall: Recall score.
            - model_f1: F1 score.
    """
    model_accuracy = accuracy_score(y_true, y_pred)
    model_precision = precision_score(y_true, y_pred, average=avg)
    model_recall = recall_score(y_true, y_pred, average=avg)
    model_f1 = f1_score(y_true, y_pred, average=avg)

    r = [model_accuracy, model_precision, model_recall, model_f1]

    return r