#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:37:53 2022

@author: elia
"""

def get_misclassified(y_true, y_pred):
    """
    Calculates the number of misclassified samples by comparing the true labels
    with the predicted labels.

    Args:
        y_true (list): List of true labels.
        y_pred (list): List of predicted labels.

    Returns:
        int: Number of misclassified samples.
    """
    return sum(1 for true, pred in zip(y_true, y_pred) if true != pred)


def misclassified_gld_nsd_sd(y_pred_golden_train, y_pred_golden_test, y_pred_nsd, y_pred_sd,
                             y_true_golden_train, y_true_golden_test, y_true_nsd, y_true_sd):
    """
    Calculates the number of misclassified samples for different datasets.

    Args:
        y_pred_golden_train (list): Predicted labels for the golden train dataset.
        y_pred_golden_test (list): Predicted labels for the golden test dataset.
        y_pred_nsd (list): Predicted labels for the NoSegDup dataset.
        y_pred_sd (list): Predicted labels for the SegDup dataset.
        y_true_golden_train (list): True labels for the golden train dataset.
        y_true_golden_test (list): True labels for the golden test dataset.
        y_true_nsd (list): True labels for the NoSegDup dataset.
        y_true_sd (list): True labels for the SegDup dataset.

    Returns:
        tuple: A tuple containing the number of misclassified samples for each dataset
        in the following order:
        - Number of misclassified samples in the golden train dataset.
        - Number of misclassified samples in the golden test dataset.
        - Number of misclassified samples in the NoSegDup dataset.
        - Number of misclassified samples in the SegDup dataset.
    """
    
    mis_gold_train = get_misclassified(y_true_golden_train, y_pred_golden_train)
    mis_gold_test = get_misclassified(y_true_golden_test, y_pred_golden_test)
    mis_nsd = get_misclassified(y_true_nsd, y_pred_nsd)
    mis_sd = get_misclassified(y_true_sd, y_pred_sd)
    
    return mis_gold_train, mis_gold_test, mis_nsd, mis_sd