#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 19:23:18 2022

@author: elia
"""

#function to plot the scores resulting from the model selection procedure
import os
import matplotlib
import matplotlib.pyplot as plt

font = {'fontname':'Arial'}
plt.rcParams.update({'font.size': 15})
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')

def plot_training_scores(cv_means, cv_stds, modelname, optimize_for, savepath, color="blue", lw=2):
    """
    Plots the training scores with standard deviations for a specific model.

    Parameters:
    cv_means (numpy.ndarray): The mean scores for each iteration.
    cv_stds (numpy.ndarray): The standard deviations of the scores for each iteration.
    modelname (str): The name of the model.
    optimize_for (str): The metric to optimize for.
    savepath (str): The path to save the generated plot.
    color (str): The color of the plot (default is "blue").
    lw (int): The line width of the plot (default is 2).

    Returns:
    The function does not return any value explicitly, but it saves the plot with the specified modelname and optimize_for in the specified savepath.
    """

    figname = "{}_training_scores.png".format(modelname)
    plt.figure()
    plt.ylim([0,1])
    plt.errorbar(range(len(cv_means)), cv_means, yerr=cv_stds, color=color, lw=lw, label="{} score".format(modelname))
    plt.xlabel('Iterations')
    plt.ylabel('{}'.format(optimize_for))
    plt.title('{} during {} training with std'.format(optimize_for, modelname), **font)
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(savepath, figname), dpi=300)
    plt.close()

    return
