#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:22:26 2022

@author: elia
"""

#function that tests the XLR autoencoder on XLR test, NSD and SD dataset
#prints performances to screen and saves to file
import os
import numpy as np

from ML_resources.get_reconstr_error import get_reconstr_error
from ML_resources.get_scaler_name import get_scaler_name
from ML_resources.plot_reconstruction_error_hist import plot_reconstruction_error_hist

def test_autoencoder(current_output_folder, ae_scaler, aenc, x_train_or, x_test_or, x_nsd_or, x_sd_or):
    """
    Tests the autoencoder by performing the following steps:
    1. Scales the input data (x_train_or, x_test_or, x_nsd_or, x_sd_or) using ae_scaler.
    2. Reconstructs the scaled test set, NSD set, and SD set using the autoencoder (aenc).
    3. Calculates the reconstruction errors for the scaled train set, test set, NSD set, and SD set.
    4. Calculates the mean, variance, and standard deviation of the reconstruction errors for each set.
    5. Plots the histograms of the reconstruction errors.
    6. Prints the mean, minimum, maximum, variance, and standard deviation of the reconstruction errors for each set.
    7. Saves the reconstruction error statistics to a log txt file.
    
    Parameters:
    - current_output_folder (str): The path to the output folder where the log txt file will be saved.
    - ae_scaler: The scaler used to scale the input data.
    - aenc: The autoencoder model used for reconstruction.
    - x_train_or: The original unscaled train set.
    - x_test_or: The original unscaled test set.
    - x_nsd_or: The original unscaled NSD set.
    - x_sd_or: The original unscaled SD set.
    
    Returns:
    - xlr_train_mean_rec_error (float): The mean reconstruction error on the scaled train set.
    - xlr_train_rec_error_std (float): The standard deviation of the reconstruction error on the scaled train set.
    """
    
    print("Scaling XLR train set, XLR test set, NSD set and SD set.")
    x_train_or_scaled = ae_scaler.transform(x_train_or)
    x_test_or_scaled = ae_scaler.transform(x_test_or)
    x_nsd_or_scaled = ae_scaler.transform(x_nsd_or)
    x_sd_or_scaled = ae_scaler.transform(x_sd_or)
    
    print("Reconstructing XLR test set, NSD set and SD set with the AE.")
    x_train_ae_scaled_reconst = aenc.predict(x_train_or_scaled)
    x_test_ae_scaled_reconst = aenc.predict(x_test_or_scaled)
    x_nsd_ae_scaled_reconst = aenc.predict(x_nsd_or_scaled)
    x_sd_ae_scaled_reconst = aenc.predict(x_sd_or_scaled)
    
    #get reconstruction errors
    print("Calculating reconstruction errors")
    x_train_rec_errors = get_reconstr_error(x_train_or_scaled, x_train_ae_scaled_reconst)    
    x_test_rec_errors = get_reconstr_error(x_test_or_scaled, x_test_ae_scaled_reconst)
    x_nsd_rec_errors = get_reconstr_error(x_nsd_or_scaled, x_nsd_ae_scaled_reconst)
    x_sd_rec_errors = get_reconstr_error(x_sd_or_scaled, x_sd_ae_scaled_reconst)
    
    print("Calculating reconstruction error mean")
    xlr_train_mean_rec_error = np.mean(x_train_rec_errors)
    xlr_test_mean_rec_error = np.mean(x_test_rec_errors)
    nsd_test_mean_rec_error = np.mean(x_nsd_rec_errors)
    sd_test_mean_rec_error = np.mean(x_sd_rec_errors)
    
    print("Calculating reconstruction error variance")
    xlr_train_rec_error_var = np.var(x_train_rec_errors)
    xlr_test_rec_error_var = np.var(x_test_rec_errors)
    nsd_test_rec_error_var = np.var(x_nsd_rec_errors)
    sd_test_rec_error_var = np.var(x_sd_rec_errors)
    
    print("Calculating reconstruction error standard deviation")
    xlr_train_rec_error_std = np.std(x_train_rec_errors)
    xlr_test_rec_error_std = np.std(x_test_rec_errors)
    nsd_test_rec_error_std = np.std(x_nsd_rec_errors)
    sd_test_rec_error_std = np.std(x_sd_rec_errors)
    
    print("Plotting reconstruction errors")
    
    plot_reconstruction_error_hist(current_output_folder, x_train_rec_errors, "XLR_train")
    plot_reconstruction_error_hist(current_output_folder, x_test_rec_errors, "XLR_test")
    plot_reconstruction_error_hist(current_output_folder, x_nsd_rec_errors, "NSD")
    plot_reconstruction_error_hist(current_output_folder, x_sd_rec_errors, "SD")
    plot_reconstruction_error_hist(current_output_folder, x_sd_rec_errors, "SD", zoom = True)
    
    xlr_train_string = "XLR train mean reconstruction error: {} | Min: {} | Max: {} | Var: {} | Std: {}".format(round(xlr_train_mean_rec_error, 2),
                                                                                                   round(min(x_train_rec_errors), 2),
                                                                                                   round(max(x_train_rec_errors), 2),
                                                                                                   round(xlr_train_rec_error_var, 2),
                                                                                                   round(xlr_train_rec_error_std, 2)
                                                                                                   )
    xlr_test_string = "XLR test mean reconstruction error: {} | Min: {} | Max: {} | Var: {} | Std: {}".format(round(xlr_test_mean_rec_error, 2),
                                                                              round(min(x_test_rec_errors), 2),
                                                                              round(max(x_test_rec_errors), 2),
                                                                              round(xlr_test_rec_error_var, 2),
                                                                              round(xlr_test_rec_error_std, 2)
                                                                              )
    
    nsd_string = "NSD mean reconstruction error: {} | Min: {} | Max: {} | Var: {} | Std: {}".format(round(nsd_test_mean_rec_error, 2),
                                                                         round(min(x_nsd_rec_errors), 2),
                                                                         round(max(x_nsd_rec_errors), 2),
                                                                         round(nsd_test_rec_error_var, 2),
                                                                         round(nsd_test_rec_error_std,2)
                                                                         ) 
    sd_string = "SD mean reconstruction error: {} | Min: {} | Max: {} | Var: {} | Std: {}".format(round(sd_test_mean_rec_error, 2),
                                                                        round(min(x_sd_rec_errors), 2),
                                                                        round(max(x_sd_rec_errors), 2),
                                                                        round(sd_test_rec_error_var, 2),
                                                                        round(sd_test_rec_error_std, 2)
                                                                        )
    
    #print to console
    print(xlr_train_string)
    print(xlr_test_string)
    print(nsd_string)
    print(sd_string)
    
    #save to log txt file
    txtname = os.path.join(current_output_folder, "XLR_{}_Autoencoder_Reconstruction_errors.txt".format(get_scaler_name("minmax")))
    f = open(os.path.join(current_output_folder, txtname), "w+")
    f.write(xlr_train_string + "\n")
    f.write(xlr_test_string + "\n")
    f.write(nsd_string + "\n")
    f.write(sd_string + "\n")
    
    #returns mean and standard deviation of reconstruction error on XLR train set for later usage
    return xlr_train_mean_rec_error, xlr_train_rec_error_std