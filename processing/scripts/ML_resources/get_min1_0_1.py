#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:35:31 2022

@author: elia
"""

#function that returns the number of -1, 0, 1 in a list
#useful to count labels in y_true for testing purposes

import numpy as np

def get_min1_0_1(y_test):
    """
    Count the occurrences of -2, -1, 0, 1, and 2 in the given NumPy array.

    Args:
        y_test (numpy.ndarray): The NumPy array containing the labels.

    Returns:
        tuple: A tuple containing the counts of -2, -1, 0, 1, and 2, respectively.
    """
    num_min2 = np.count_nonzero(y_test == -2)
    num_min1 = np.count_nonzero(y_test == -1)
    num_0 = np.count_nonzero(y_test == 0)
    num_1 = np.count_nonzero(y_test == 1)
    num_2 = np.count_nonzero(y_test == 2)

    return num_min2, num_min1, num_0, num_1, num_2

