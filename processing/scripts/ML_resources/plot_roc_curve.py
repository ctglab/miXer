#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 23:19:16 2022

@author: elia
"""

#function that plots the precision recall curve
#can handle both 2 class problems and 3 class problems
#3 (4) class problems are treated as three (four) separate 1 vs all problems
import os
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

from ML_resources.onevsall_conv import onevsall_conv

font = {'fontname':'Arial'}
plt.rcParams.update({'font.size': 15})
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')

def plot_roc_curve(modelname, savepath, y_true, y_pred):
    """
    Plots the Receiver Operating Characteristic (ROC) curve for a classification model.

    Parameters:
    modelname (str): The name of the model.
    savepath (str): The path to save the generated plot.
    y_true (numpy.ndarray): The true labels.
    y_pred (numpy.ndarray): The predicted probabilities.

    Returns:
    The function does not return any value explicitly, but it saves the plot with the specified modelname in the specified savepath.
    """
    
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    lw = 2
    
    classes = list(set(y_true))
    # Original RGB color values
    ddel_col = (0, 137, 255)
    del_col = (0, 213, 255)
    wt_col = (128, 128, 128)
    dup_col = (255, 128, 0)
    mdup_col = (255, 0, 0)
    
    # Convert to 0-1 range
    ddel_col = (ddel_col[0] / 255, ddel_col[1] / 255, ddel_col[2] / 255)
    del_col = (del_col[0] / 255, del_col[1] / 255, del_col[2] / 255)
    wt_col = (wt_col[0] / 255, wt_col[1] / 255, wt_col[2] / 255)
    dup_col = (dup_col[0] / 255, dup_col[1] / 255, dup_col[2] / 255)
    mdup_col = (mdup_col[0] / 255, mdup_col[1] / 255, mdup_col[2] / 255)
    
    #Adapting for XLR vs NSD/SD dataset classes
    if len(classes) > 3:
        color_list =[ddel_col, del_col, wt_col, dup_col, mdup_col]
    else:
        color_list = [del_col, wt_col, dup_col]

    #TODO non è più un problema a due classi, sistemare
    if len(classes) == 2:
        #stampa solo la classe positiva
        fpr, tpr, _ = roc_curve(y_true, y_pred[:, 1], pos_label= 1)
        roc_auc = auc(fpr, tpr)
        
        plt.figure()
        plt.plot(fpr, tpr, color = "red", lw = lw, label = "ROC curve class 1 (AUC = {:.2f})".format(roc_auc))
        plt.plot([0,1], [0,1], color = "gray", lw = lw, linestyle = "--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', **font)  
        plt.ylabel('True Positive Rate', **font) 
        plt.title('ROC {}'.format(modelname), **font)  
        plt.legend(loc="lower right")

    else:
        #concateno y_true, y_pred
        y_reshape = y_true.reshape(-1,1)
        y_stack = np.hstack((y_reshape,y_pred))
        
        colors = itertools.cycle(color_list)
        plt.figure()
        for c, color in zip(classes, colors):
            #le righe dove la predizione è uguale alla classe corrente c vengono messe uguale ad 1, le altre a 0
            #la probabilità viene lasciata uguale dove la predizione è uguale a c
            #altimenti viene sostituita con (1 - il suo valore) (probabilità dell'evento complementare)
            y_converted = onevsall_conv(y_stack, c, classes)
            #funziona, ora plotto la roc per ogni classe come se fosse la classe positiva!
            fpr, tpr, _ = roc_curve(y_converted[:,0].astype("int"), y_converted[:,1], pos_label = 1)
            roc_auc = auc(fpr, tpr)
            plt.plot(fpr, tpr, color = color, lw = lw, label = "ROC class {} (AUC = {:.2f})".format(c,roc_auc))
        
        plt.plot([0,1], [0,1], color = "gray", lw = lw, linestyle = "--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', **font)  
        plt.ylabel('True Positive Rate', **font) 
        plt.title('ROC {}'.format(modelname), **font)  
        plt.legend(loc="lower right")
        
    figname = "{}_ROC.pdf".format(modelname)
    plt.savefig(os.path.join(savepath, figname), dpi = 300)
    plt.close()
    return