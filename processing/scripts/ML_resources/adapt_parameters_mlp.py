#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:21:10 2022

@author: elia
"""

#function to adapt the input parameters to scikit learn's wishes
from ML_resources.str_to_bool import str_to_bool

def adapt_parameters_mlp(params):
    """
    Adapt the parameters for training a multilayer perceptron (MLP) model.

    The function modifies the input `params` dictionary to adapt the parameters
    for training an MLP model. It generates the appropriate `hidden_layer_sizes`
    parameter based on the provided `hidden_layer_units` and `hidden_layers` lists.
    It also handles the conversion of the `early_stopping` parameter from a boolean
    string representation to a boolean value.

    Args:
        params (dict): A dictionary containing the parameters for training the MLP model.
            The dictionary should include the following keys:
                - "hidden_layer_units": A list of integers representing the number of units
                  in each hidden layer.
                - "hidden_layers": A list of integers representing the number of hidden layers.
                - "early_stopping": A list containing a boolean string representation (e.g., "True" or "False")
                  indicating whether early stopping should be used.

    Returns:
        dict: The modified `params` dictionary with the following updates:
            - "hidden_layer_sizes": A list of tuples, where each tuple represents the number of units
              in each hidden layer, based on the combinations of "hidden_layer_units" and "hidden_layers".
            - "early_stopping": A list containing a boolean value indicating whether early stopping should be used.
    """
    units_list = params["hidden_layer_units"]
    hlayers_list = params["hidden_layers"]
    early_stopping = params["early_stopping"][0] #TODO non riesco a passare un booleano con JSON, workaround terribile per ora
    
    del params["early_stopping"]
    
    if str_to_bool(early_stopping, "Early stopping") == True:
        params["early_stopping"] = [True]
    else:
        params["early_stopping"] = [False]
        
    #remove from dictionary
    del params["hidden_layer_units"] #this mutates the existing dictionary, so the contents of the dictionary changes for anybody else who has a reference to the same instance
    del params["hidden_layers"]
    
    #build the tuples
    tuples_list = []
    tmp = []
    
    for units in units_list:
        for layers in hlayers_list:
            for i in range(layers):
                tmp.append(units)
            tuples_list.append(tuple(tmp))
            tmp = []
    
    params["hidden_layer_sizes"] = tuples_list
    return params
