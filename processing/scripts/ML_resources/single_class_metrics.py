#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 23:50:25 2024

@author: elia
"""

import pandas as pd
import numpy as np
from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix

def single_class_metrics(df, classes, modelname):
    # Extracting true labels and predictions
    lbl = "{}_pred".format(modelname.upper())
    true_labels = df['Class']
    predicted_labels = df[lbl]
    
    # List to store metrics for each class
    precision_per_class = []
    recall_per_class = []
    f1_per_class = []
    tpr_per_class = []
    tnr_per_class = []
    fpr_per_class = []
    fnr_per_class = []
    
    # Loop over each class
    for c in classes:
        # Converting the problem to binary classification
        true_labels_binary = (true_labels == c).astype(int)
        predicted_labels_binary = (predicted_labels == c).astype(int)
        
        # Calculating metrics for the current class
        precision_c = precision_score(true_labels_binary, predicted_labels_binary)
        recall_c = recall_score(true_labels_binary, predicted_labels_binary)
        f1_c = f1_score(true_labels_binary, predicted_labels_binary)
        
        # Confusion matrix
        tn, fp, fn, tp = confusion_matrix(true_labels_binary, predicted_labels_binary).ravel()
        
        # Calculating TPR, TNR, FPR, FNR
        tpr = tp / (tp + fn) if (tp + fn) != 0 else 0
        tnr = tn / (tn + fp) if (tn + fp) != 0 else 0
        fpr = fp / (fp + tn) if (fp + tn) != 0 else 0
        fnr = fn / (fn + tp) if (fn + tp) != 0 else 0
        
        # Appending metrics to the respective lists
        precision_per_class.append(precision_c)
        recall_per_class.append(recall_c)
        f1_per_class.append(f1_c)
        tpr_per_class.append(tpr)
        tnr_per_class.append(tnr)
        fpr_per_class.append(fpr)
        fnr_per_class.append(fnr)
    
    # Create a DataFrame
    metrics_df = pd.DataFrame({
        "Precision": precision_per_class,
        "Recall": recall_per_class,
        "F1 Score": f1_per_class,
        "True Positive Rate": tpr_per_class,
        "True Negative Rate": tnr_per_class,
        "False Positive Rate": fpr_per_class,
        "False Negative Rate": fnr_per_class
    }, index=classes)
    
    return metrics_df