#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 09:38:41 2022

@author: elia
"""

#function that plots the histogram of prediction probability by class
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')
plt.rcParams.update({'font.size': 15})

def plot_pred_proba(modelname, testname, dataset, savepath, bins = 100, dpi = 300):
    """
    Plots the histogram of predicted probabilities for a model.

    Args:
        modelname (str): Name of the model.
        testname (str): Name of the test.
        dataset (pandas.DataFrame): Dataset containing the predicted probabilities.
        savepath (str): Path to save the generated plot.
        bins (int, optional): Number of bins in the histogram. Defaults to 100.
        dpi (int, optional): Resolution of the saved plot in dots per inch. Defaults to 300.

    Returns:
        None
    """
    
    fig, ax = plt.subplots(figsize = (8,6))
    pred_column_name = "{}_pred_proba_".format(modelname)
    heights_min1, bins_min1 = np.histogram(dataset[pred_column_name + "(-1)"])
    heights_0, bins_0 = np.histogram(dataset[pred_column_name + "(0)"], bins = bins_min1)
    heights_1, bins_1 = np.histogram(dataset[pred_column_name + "(1)"], bins = bins_min1)
    
    width = (bins_min1[1] - bins_min1[0])/5
    
    bar_m1 = ax.bar(bins_min1[:-1] - (width), heights_min1, width = width, facecolor = "blue")
    bar_0 = ax.bar(bins_0[:-1], heights_0, width = width, facecolor = "orange")
    bar_1 = ax.bar(bins_1[:-1]+ (width), heights_1, width = width, facecolor = "green")
    
    plt.xlabel("Model confidence")
    plt.ylabel("Occurrencies (blue = -1, orange = 0, green = 1)")
    plt.title("{}_{}_proba_hist".format(testname, modelname))
    figtitle = "{}_{}_proba_hist.png".format(testname, modelname)
    
    fig.savefig(os.path.join(savepath, figtitle), dpi = dpi)
    plt.close()

    return