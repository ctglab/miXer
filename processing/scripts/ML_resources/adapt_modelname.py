#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:18:15 2022

@author: elia
"""

#function that adapts the script's model name for saving purposes

def adapt_modelname(modelname):
    """
    Adapt the model name to a standardized format.
    
    Parameters:
    - modelname (str): The model name or identifier.
    
    Returns:
    - adapted_modelname (str): The standardized model name.
    
    Raises:
    - ValueError: If the provided model name is not supported.
    """
    
    if modelname.lower() in ["svc", "svm"]:
        return "SVC"
    elif modelname.lower() in ["rf", "random forest", "randomforest"]:
        return "RF"
    else:
        raise ValueError("Unsupported model identifier") #TODO test
        
    