#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:10:45 2022

@author: elia
"""

#function that given a scaler identifier, returns the correct scikit learn object
from sklearn.preprocessing import  RobustScaler, StandardScaler, MinMaxScaler
import warnings

def get_scaler(scaler_id, verbose = 0):
    """
    Get a scaler object based on the scaler identifier.

    Args:
        scaler_id (str): Scaler identifier.
        verbose (int, optional): Verbosity level. Set to 0 for no output. Defaults to 0.

    Returns:
        object: Scaler object.

    Raises:
        Warning: If an unavailable scaler identifier is provided, it defaults to RobustScaler.

    Examples:
        >>> scaler_id = "rs"
        >>> scaler = get_scaler(scaler_id)
    """
    
    if scaler_id.lower() in ["robustscaler", "rs", "robust", "robust_scaler"]:
        if verbose > 0:
            print("Using Robust Scaler")
        return RobustScaler()
    elif scaler_id.lower() in ["standardscaler", "ss", "stdsc", "standard_scaler"]:
        if verbose > 0:
            print("Using Standard Scaler")
        return StandardScaler()
    elif scaler_id.lower() in ["minmaxscaler", "minmax", "mmx", "minmax_scaler", "mmscaler"]:
        if verbose > 0:
            print("Using MinMax Scaler")
        return MinMaxScaler()
    else:
        warnings.warn("Unavailable scaler identifier, defaulting to RobustScaler", Warning) #Should not happen since get_scaler_name.py outputs a correct identifier. TODO error and halt exec.
        return RobustScaler()