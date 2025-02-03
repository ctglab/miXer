#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 23:21:59 2022

@author: elia
"""

#function that converts a multiclass problem in one-vs-all problems
import numpy as np


def onevsall_conv(y_stack, c, classes):
    """
    Converts the multi-class labels to one-vs-all format and computes the predicted probabilities for a specific class.

    Parameters:
    y_stack (numpy.ndarray): The input multi-class labels stacked in a 2D array.
    c (int): The target class for which the conversion and probabilities are computed.
    classes (list): The list of class labels.

    Returns:
    numpy.ndarray: The converted labels and predicted probabilities for the specified target class.
    """

    rows = y_stack.shape[0]
    out = np.zeros((y_stack.shape[0], 2))  # column for y_true, y_pred_proba

    for i in range(rows):
        y_row = y_stack[i]
        pred_index = np.argmax(y_row[1:])  # indicates the model's output. E.g., pred_index = 1 means the model's highest confidence is for class -1 (del), pred_index = 0 for class -2, pred_index = 2 for class 0, pred_index = 3 for class 1 (dup)

        if int(y_row[0]) == c:
            out[i][0] = 1  # considered as a positive class, so assign 1
        else:
            out[i][0] = 0  # can be any of the other two classes, so considered negative, 0 (one-vs-all classifier)

        if c == (pred_index + min(classes)):  # an ugly hack, but ran out of ideas. We assign the classes as [-2, -1, 0, 1, 2], the index of prediction can be [0, 1, 2, 3,4], so simply subtract 2 from the index
            out[i][1] = y_row[pred_index + 1]  # assigned probability is the one that the model assigns to the positive class
        else:  # +1 is the offset due to the fact that the first column represents y_true
            out[i][1] = 1 - y_row[pred_index + 1]  # assigned probability is the sum of the probabilities that the model assigns to the negative classes

    return out