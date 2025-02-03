#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 18:20:49 2022

@author: elia
"""

#function to plot AE reconstruction error histograms
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')
plt.rcParams.update({'font.size': 15})

def plot_reconstruction_error_hist(current_output_folder, rec_err, test_type, add_mean_std = True, zoom = False, n_bins = None):
    """
    Plots a histogram of the reconstruction errors for an autoencoder.

    Args:
        current_output_folder (str): Path to the current output folder.
        rec_err (numpy.ndarray): Array of reconstruction errors.
        test_type (str): Type of the test.
        add_mean_std (bool, optional): Whether to add mean and standard deviation lines to the plot. Defaults to True.
        zoom (bool, optional): Whether to zoom in on the x-axis. Defaults to False.
        n_bins (int, optional): Number of bins in the histogram. If None, the Freedman–Diaconis rule is used to select the number of bins. Defaults to None.

    Returns:
        None
    """
    
    #if no specific number of bins is specified
    #this will follow the Freedman–Diaconis rule to select the number of bins in an histogram
    if n_bins is None:
        q25, q75 = np.percentile(rec_err, [25, 75])
        bin_width = 2 * (q75 - q25) * len(rec_err) ** (-1/3)
        n_bins = int(round((max(rec_err) - min(rec_err)) / bin_width))
        #print("Freedman–Diaconis number of bins:", bins)
        
    plt.hist(rec_err, density=False, bins=n_bins)
    
    plt.ylabel("Target regions")
    plt.xlabel("AE Reconstruction Error")
    
    if add_mean_std:
        mean_re = np.mean(rec_err)
        std_re = np.std(rec_err)
        
        plt.axvline(x = mean_re, label='Mean RE = {}'.format(round(mean_re,2)), c= "red")
        plt.axvline(x = mean_re + std_re, label='Mean RE + STD = {}'.format(round(mean_re + std_re,2)), c= "green")
        #plt.axvline(x = mean_re - std_re, label='Mean RE - STD = {}'.format(round(mean_re - std_re,2)), c= "green")
        plt.axvline(x = mean_re + 2*std_re, label='Mean RE + 2xSTD = {}'.format(round(mean_re + 2*std_re,2)), c= "orange")
        #plt.axvline(x = mean_re - 2*std_re, label='Mean RE - 2xSTD = {}'.format(round(mean_re - 2*std_re,2)), c= "orange")
        
    if zoom:
        plt.xlim(0, 5)
        figname = "AE_{}_reconstruction_error_hist_zoom.png".format(test_type)
        
    else:
        figname = "AE_{}_reconstruction_error_hist.png".format(test_type)
    
    plt.legend()
    plt.savefig(os.path.join(current_output_folder, figname), dpi = 300)
    plt.close()
    
    return