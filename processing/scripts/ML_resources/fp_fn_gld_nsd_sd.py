#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:45:05 2022

@author: elia
"""

import numpy as np
#functions that calculate the number of TP and FN 

def fp_fn(y_true, y_pred):
    """
    Calculate the number of false positives and false negatives between the true labels and predicted labels,
    along with True Positive Rate (TPR), False Positive Rate (FPR), True Negative Rate (TNR), and False Negative Rate (FNR).

    Args:
        y_true (array-like): The true labels.
        y_pred (array-like): The predicted labels.

    Returns:
        tuple: A tuple containing the number of false positives, false negatives, TPR, FPR, TNR, and FNR, respectively.
    """
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    # Find false positives and false negatives using boolean indexing
    fp = np.sum((y_true == 0) & (y_pred != 0))
    fn = np.sum((y_true != 0) & (y_pred == 0))
    
    # Calculate True Positive Rate (TPR), False Positive Rate (FPR),
    # True Negative Rate (TNR), and False Negative Rate (FNR)
    tp = np.sum((y_true != 0) & (y_pred != 0))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    tpr = tp / (tp + fn) if (tp + fn) != 0 else 0
    fpr = fp / (fp + tn) if (fp + tn) != 0 else 0
    tnr = tn / (tn + fp) if (tn + fp) != 0 else 0
    fnr = fn / (fn + tp) if (fn + tp) != 0 else 0
    
    return fp, fn, tpr, fpr, tnr, fnr


def fp_fn_gld_nsd_sd(y_pred_golden_train, y_pred_golden_test, y_pred_nsd, y_pred_sd,
                     y_true_golden_train, y_true_golden_test, y_true_nsd, y_true_sd):
    """
    Calculate the number of false positives and false negatives for different datasets.

    Args:
        y_pred_golden_train (array-like): The predicted labels for the golden train set.
        y_pred_golden_test (array-like): The predicted labels for the golden test set.
        y_pred_nsd (array-like): The predicted labels for the NSD set.
        y_pred_sd (array-like): The predicted labels for the SD set.
        y_true_golden_train (array-like): The true labels for the golden train set.
        y_true_golden_test (array-like): The true labels for the golden test set.
        y_true_nsd (array-like): The true labels for the NSD set.
        y_true_sd (array-like): The true labels for the SD set.

    Returns:
        tuple: A tuple containing the results of fp_fn function for each dataset,
            in the following order: 
            - For golden train set: [fp_gold_train, fn_gold_train, tpr_gold_train, fpr_gold_train, tnr_gold_train, fnr_gold_train]
            - For golden test set: [fp_gold_test, fn_gold_test, tpr_gold_test, fpr_gold_test, tnr_gold_test, fnr_gold_test]
            - For NSD set: [fp_nsd, fn_nsd, tpr_nsd, fpr_nsd, tnr_nsd, fnr_nsd]
            - For SD set: [fp_sd, fn_sd, tpr_sd, fpr_sd, tnr_sd, fnr_sd]
    """
    fp_gold_train, fn_gold_train, tpr_gold_train, fpr_gold_train, tnr_gold_train, fnr_gold_train = fp_fn(y_true_golden_train, y_pred_golden_train)
    fp_gold_test, fn_gold_test, tpr_gold_test, fpr_gold_test, tnr_gold_test, fnr_gold_test = fp_fn(y_true_golden_test, y_pred_golden_test)
    fp_nsd, fn_nsd, tpr_nsd, fpr_nsd, tnr_nsd, fnr_nsd = fp_fn(y_true_nsd, y_pred_nsd)
    fp_sd, fn_sd, tpr_sd, fpr_sd, tnr_sd, fnr_sd = fp_fn(y_true_sd, y_pred_sd)
    
    
    return ([fp_gold_train, fn_gold_train, tpr_gold_train, fpr_gold_train, tnr_gold_train, fnr_gold_train],
            [fp_gold_test, fn_gold_test, tpr_gold_test, fpr_gold_test, tnr_gold_test, fnr_gold_test],
            [fp_nsd, fn_nsd, tpr_nsd, fpr_nsd, tnr_nsd, fnr_nsd],
            [fp_sd, fn_sd, tpr_sd, fpr_sd, tnr_sd, fnr_sd])


