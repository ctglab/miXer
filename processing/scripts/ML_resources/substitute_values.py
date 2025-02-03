#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 23:13:47 2023

@author: elia
"""

def substitute_values(vector):
    """
    Substitutes values in a NumPy vector.
    
    Replaces -2 values with -1 and values greater or equal than 2 with +1 in the given vector.

    Parameters:
    vector (numpy.ndarray): The input vector to be modified.

    Returns:
    numpy.ndarray: The modified vector with substituted values.
    """
    
    vector[vector == -2] = -1
    vector[vector >= 2] = 1
    return vector