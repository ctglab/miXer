#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 23:14:05 2022

@author: elia
"""

#function that plots a confusion matrix
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import itertools
import os

font = {'fontname':'Arial'}
plt.rcParams.update({'font.size': 15})
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')

def plot_confusion_matrix(modelname, cm, classes, savepath, cmap=plt.cm.Blues):
    """
    Prints and plots the confusion matrix for a given model.

    Args:
        modelname (str): Name of the model.
        cm (numpy.ndarray): Confusion matrix.
        classes (list): List of class labels.
        savepath (str): Path to save the generated plot.
        cmap (matplotlib colormap, optional): Colormap for the plot. Defaults to plt.cm.Blues.

    Returns:
        None
    """
    figtitle = "{} Confusion Matrix".format(modelname)
    figname = "{}_Confusion_matrix.pdf".format(modelname)
   
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(figtitle)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes)
    plt.yticks(tick_marks, classes)

    fmt = 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center",
                 color="white" if cm[i, j] > thresh else "black",
                 **font)

    plt.ylabel('True label', **font)
    plt.xlabel('Predicted label', **font)
    plt.tight_layout()
    plt.savefig(os.path.join(savepath, figname), dpi=300)
    plt.close()
    return
