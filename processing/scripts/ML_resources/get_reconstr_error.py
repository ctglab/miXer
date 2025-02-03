#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:27:48 2022

@author: elia
"""

#function that returns a vector containing the reconstruction error
#reconstruction error is calculated as the Euclidean norm between true and reconstructed vectors
import numpy as np

def get_reconstr_error(x_true, x_reconstr):
    """
    Calculate the reconstruction error between true samples and reconstructed samples.

    Args:
        x_true (numpy.ndarray): Array of true samples.
        x_reconstr (numpy.ndarray): Array of reconstructed samples.

    Returns:
        list: List of reconstruction errors for each sample.

    Examples:
        >>> true_samples = np.array([[1, 2, 3], [4, 5, 6]])
        >>> reconstructed_samples = np.array([[0.9, 2.1, 3.2], [3.9, 5.2, 6.1]])
        >>> errors = get_reconstr_error(true_samples, reconstructed_samples)
    """
    
    errors = []
    for i in range(len(x_true)):
        t = x_true[i]
        r = x_reconstr[i]
        errors.append(np.linalg.norm(t-r))
    return errors