#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:21:26 2022

@author: elia
"""

#function to plot an MLP
#uses VisualizeNN from Jianzheng Liu (further info and author contact info in VisualizeNN script) to visualize the MLP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import ML_resources.VisualizeNN as VisNN

#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')
plt.rcParams.update({'font.size': 15})

def plot_nn(x_train_scaled, nn, out_folder, disp_columns, modelname, nn_accuracy):
    """
    Plots the neural network architecture of the best found MLP model.

    Args:
        x_train_scaled (numpy.ndarray): Scaled input training data.
        nn (GridSearchCV): GridSearchCV object containing the trained MLP model.
        out_folder (str): Output folder path for saving the network visualization.
        disp_columns (int): Number of columns to display in the network visualization.
        modelname (str): Name of the model.
        nn_accuracy (float): Accuracy of the MLP model.

    Returns:
        None
    """
    
    #plot the mlp if needed
    print("Drawing the best found MLP")
    print("Output layer activation: {} | Number of units in output layer: {}".format(nn.best_estimator_.out_activation_, nn.best_estimator_.n_outputs_))
    network_structure = np.hstack(([x_train_scaled.shape[1]], np.asarray(nn.best_params_["hidden_layer_sizes"]), nn.best_estimator_.n_outputs_))

    # Draw the Neural Network with weights
    network=VisNN.DrawNN(network_structure, nn.best_estimator_.coefs_, out_folder, modelname, nn_accuracy, 0, disp_columns)
    network.draw()
    plt.close()
    return
