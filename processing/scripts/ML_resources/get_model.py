#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:19:45 2022

@author: elia
"""

#function that returns the correct scikit learn object given a model identifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

def get_model(identifier):
    """
    Get an instance of the specified model.

    Args:
        identifier (str): The identifier of the model.

    Returns:
        object: An instance of the specified model.

    Raises:
        ValueError: If the model identifier is not supported.

    Examples:
        >>> model = get_model("svc")
    """
    
    if identifier.lower() in ["svc", "svm"]:
        return SVC()
    elif identifier.lower() in ["rf", "rfc"]:
        return RandomForestClassifier()
    else:
        raise ValueError("Unsupported model identifier.")