#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:54:31 2022

@author: elia
"""

#function to plot the autoencoder
#uses VisualizeNN from Jianzheng Liu (further info and author contact info in VisualizeNN script) to visualize the MLP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import ML_resources.VisualizeNN as VisNN

#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')

def plot_ae(x_train, ae, out_folder, disp_columns, ae_name, reconstr_error):
    """
    Plots the Autoencoder network with weights.

    Args:
        x_train (numpy.ndarray): Input training data.
        ae (Autoencoder): Trained Autoencoder model.
        out_folder (str): Output folder path for saving the network visualization.
        disp_columns (int): Number of columns to display in the network visualization.
        ae_name (str): Name of the Autoencoder.
        reconstr_error (float): Reconstruction error.

    Returns:
        None
    """
    
    print("Drawing the Autoencoder")
    network_structure = np.hstack(([x_train.shape[1]], np.asarray(ae.hidden_layer_sizes),[x_train.shape[1]]))
    # Draw the Neural Network with weights
    network=VisNN.DrawNN(network_structure, ae.coefs_, out_folder, ae_name, 0, reconstr_error, disp_columns)
    network.draw()
    plt.close()
    return