#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:15:52 2022

@author: elia
"""

#function that returns a scikit learn autoencoder object with the requested architecture
from sklearn.neural_network import MLPRegressor

import os
import random
import numpy as np

#controlling randomness
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value) #fixed ambient variable value
random.seed(seed_value) #fixed random module seed
np.random.seed(seed_value) #fixed numpy's random module seed


def get_autoencoder(seed_value, input_shape, hidden_layer_shape=3, training_iter=10):
    """
    Create an autoencoder model using a multi-layer perceptron regressor.

    Args:
        seed_value (int): The seed value for random number generation.
        input_shape (int): The shape of the input layer.
        hidden_layer_shape (int, optional): The shape of the hidden layer. Defaults to 3.
        training_iter (int, optional): The number of training iterations. Defaults to 10.

    Returns:
        MLPRegressor: An instance of MLPRegressor representing the autoencoder model.
    """
    mlp = MLPRegressor(hidden_layer_sizes=(hidden_layer_shape + 1,
                                           hidden_layer_shape,
                                           hidden_layer_shape + 1),
                       activation="relu",
                       solver="adam",
                       learning_rate_init=0.0001,
                       max_iter=training_iter,
                       tol=0.0001,
                       random_state=seed_value,
                       verbose=True)
    return mlp
