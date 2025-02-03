#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:09:30 2022

@author: elia
"""

#function to get formatted scaler name
import warnings

def get_scaler_name(scaler_id):
    """
    Get the full name of the scaler based on the scaler identifier.

    Args:
        scaler_id (str): Scaler identifier.

    Returns:
        str: Full name of the scaler.

    Raises:
        Warning: If an unavailable scaler identifier is provided, it defaults to RobustScaler.

    Examples:
        >>> scaler_id = "rs"
        >>> scaler_name = get_scaler_name(scaler_id)
    """
    
    if scaler_id.lower() in ["robustscaler", "rs", "robust", "robust_scaler"]:
        return "robustscaler"
    elif scaler_id.lower() in ["standardscaler", "ss", "stdsc", "standard_scaler"]:
        return "standardscaler"
    elif scaler_id.lower() in ["minmaxscaler", "minmax", "mmx", "minmax_scaler", "mmscaler"]:
        return "minmaxscaler"
    else:
        warnings.warn("Unavailable scaler identifier, defaulting to RobustScaler", Warning)
        return "robustscaler"