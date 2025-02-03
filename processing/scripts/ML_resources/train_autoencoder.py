#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 15:00:21 2022

@author: elia
"""

#train or load the autoencoder for anomaly detection
import os
from joblib import load, dump

from ML_resources.get_scaler_name import get_scaler_name
from ML_resources.get_scaler import get_scaler
from ML_resources.get_autoencoder import get_autoencoder
from ML_resources.prepare_training import add_noise_to_data

def train_autoencoder(current_output_folder, x_train_or, retrain_ae, seed_value, denoising = False, training_iter = 10):
    """
    Trains an autoencoder model and saves the scaler and the autoencoder model.

    Parameters:
    - current_output_folder (str): The path to the output folder where the scaler and autoencoder model will be saved.
    - x_train_or (array-like): The input data for training the autoencoder.
    - retrain_ae (bool): Indicates whether to force retraining of the autoencoder.
    - seed_value (int): The seed value for reproducible results.
    - denoising (bool): Indicates whether the autoencoder is trained for denoising.
    - training_iter (int): The number of training iterations for the autoencoder.

    Returns:
    - ae_scaler (object): The scaler object used for scaling the data.
    - aenc (object): The trained autoencoder model.
    - ae_output_folder (str): The path to the output folder where the scaler and autoencoder model are saved.
    """
    
    #getting scaler to scale the samples
    #defining autoencoder scaler savepath
    if denoising:
        ae_folder_name = "Autoencoder_denoising"
        ae_output_folder = os.path.join(current_output_folder, ae_folder_name)
        
        ae_scalerpath = os.path.join(ae_output_folder, "AE_Denoise_" + get_scaler_name("minmax") + ".joblib")
        #defining autoencoder savepath
        aenc_savepath = os.path.join(ae_output_folder, "XLR_{}_Autoencoder_denoise.joblib".format(get_scaler_name("minmax")))
        
    else:    
        ae_folder_name = "Autoencoder_clean"
        ae_output_folder = os.path.join(current_output_folder, ae_folder_name)
        
        ae_scalerpath = os.path.join(ae_output_folder, "AE_" + get_scaler_name("minmax") + ".joblib")
        #defining autoencoder savepath
        aenc_savepath = os.path.join(ae_output_folder, "XLR_{}_Autoencoder.joblib".format(get_scaler_name("minmax")))
    
    if not os.path.exists(ae_output_folder):
        os.makedirs(ae_output_folder)
            
    if os.path.isfile(ae_scalerpath):
        print("Reloading Autoencoder's Scaler.")
        #reload scaler
        ae_scaler = load(ae_scalerpath)
    else:
        print("Training Autoencoder Scaler.")
        #train scaler
        ae_scaler = get_scaler("minmax")
        ae_scaler.fit(x_train_or)
        dump(ae_scaler, ae_scalerpath)
    
    if denoising:
        x_input = ae_scaler.transform(add_noise_to_data(0.5, 0.5, x_train_or))
        x_to_reconstruct = ae_scaler.transform(x_train_or)
    
    else:
        x_input = x_to_reconstruct = ae_scaler.transform(x_train_or)
    
    if os.path.isfile(aenc_savepath):
        if retrain_ae:
            print("Forcing Re-training of the Autoencoder.")
            aenc = get_autoencoder(seed_value, x_input.shape[1], hidden_layer_shape = 3, training_iter = training_iter)
            aenc.fit(x_input, x_to_reconstruct)
            dump(aenc, aenc_savepath)
        else:
            print("Reloading trained Autoencoder.")
            aenc = load(aenc_savepath)
    else:
        print("Training Autoencoder.")
        aenc = get_autoencoder(seed_value, x_input.shape[1], hidden_layer_shape = 3, training_iter = training_iter)
        aenc.fit(x_input, x_to_reconstruct)
        dump(aenc, aenc_savepath)
    
    return ae_scaler, aenc, ae_output_folder


