#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 13:24:44 2022

@author: elia
"""

#wrapping function for splitting datasets in samples - labels
#splits golden set in gld_training / gld_test according to specific sample number or split fraction (eg. 13k samples or 80% train-20% test)
import warnings
import os
import pandas as pd
from ML_resources.prepare_training import prepare_training

def dataset_setup(gld, nsd, sd, train_samples, split_fraction, mu, sigma, add_noise, seed_value,
                  skip_training = False, save_results = False, out_dir = None):
    """
    Prepare the dataset for training by splitting the golden dataset and optionally adding noise.
    
    Parameters:
    - gld (array-like): The golden dataset.
    - nsd (array-like): The noSegDup dataset.
    - sd (array-like): The SegDup dataset.
    - train_samples (int): The number of samples to use for training. If None, the dataset is split based on split_fraction.
    - split_fraction (float): The fraction of the dataset to use for training. Only used if train_samples is None.
    - mu (float): The mean of the noise distribution.
    - sigma (float): The standard deviation of the noise distribution.
    - add_noise (bool): Indicates whether to add noise to the dataset.
    - seed_value (int): The seed value for reproducible results.
    - skip_training (bool): Indicates whether to skip preparing the noSegDup and SegDup datasets.
    
    Returns:
    - output_dict (dict): A dictionary containing the prepared datasets and their corresponding labels.
      - "x_train_or" (array-like): The input features of the training dataset.
      - "x_test_or" (array-like): The input features of the testing dataset.
      - "y_train" (array-like): The labels of the training dataset.
      - "y_test" (array-like): The labels of the testing dataset.
      - "x_nsd_or" (array-like): The input features of the noSegDup dataset.
      - "x_sd_or" (array-like): The input features of the SegDup dataset.
      - "y_nsd" (array-like): The labels of the noSegDup dataset.
      - "y_sd" (array-like): The labels of the SegDup dataset.
    """
    
    if train_samples is None:
        output_dict = prepare_training(gld, split_to_train = True,
                                       mu = mu, sigma = sigma,
                                       split_fraction = split_fraction,
                                       add_noise = add_noise,
                                       seed_value = seed_value
                                       )
        
       
    else:
        output_dict = prepare_training(gld, split_to_train = True,
                                       mu = mu, sigma = sigma,
                                       train_samples = train_samples,
                                       add_noise = add_noise,
                                       seed_value = seed_value
                                       )
    
    if not skip_training:
        #prepare noSegDup and SegDup datasets
        x_nsd_or, y_nsd = prepare_training(nsd, seed_value = seed_value)
        x_sd_or, y_sd = prepare_training(sd, seed_value = seed_value)
        
        #adding nosegdup and segdup to output_dict
        output_dict["x_nsd_or"] = x_nsd_or
        output_dict["x_sd_or"] = x_sd_or
        output_dict["y_nsd"] = y_nsd
        output_dict["y_sd"] = y_sd
    
    #check to save the resulting datasets
    if save_results:
        print("Saving prepared ChrX datasets")
        #building dataframes
        #XLR_train dataframe -> noisy - noise only - clean only
        ds_xlr_train_clean = pd.concat([pd.DataFrame(output_dict["x_train_clean_only"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                    pd.Series(output_dict["y_train_clean_only"], name='Class')], axis=1)
        ds_xlr_test_clean = pd.concat([pd.DataFrame(output_dict["x_test_clean_only"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                    pd.Series(output_dict["y_test_clean_only"], name='Class')], axis=1)
        #saving datasets
        ds_xlr_train_clean.to_csv(os.path.join(out_dir,"XLR_train_clean_only.txt.gz"), compression = "gzip", sep = "\t", index = False)
        ds_xlr_test_clean.to_csv(os.path.join(out_dir,"XLR_test_clean_only.txt.gz"), compression = "gzip", sep = "\t", index = False)
        
        if add_noise:
            ds_xlr_train_noisy_only = pd.concat([pd.DataFrame(output_dict["x_train_noisy_only"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_train_noisy_only"], name='Class')], axis=1)
            ds_xlr_test_noisy_only = pd.concat([pd.DataFrame(output_dict["x_test_noisy_only"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_test_noisy_only"], name='Class')], axis=1)
            
            ds_xlr_train_with_noise = pd.concat([pd.DataFrame(output_dict["x_train_with_noise"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_train_with_noise"], name='Class')], axis=1)
            ds_xlr_test_with_noise = pd.concat([pd.DataFrame(output_dict["x_test_with_noise"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_test_with_noise"], name='Class')], axis=1)
            
            ds_xlr_train_noisy_only.to_csv(os.path.join(out_dir,"XLR_train_noise_only.txt.gz"), compression = "gzip", sep = "\t", index = False)
            ds_xlr_test_noisy_only.to_csv(os.path.join(out_dir,"XLR_test_noise_only.txt.gz"), compression = "gzip", sep = "\t", index = False)
            
            ds_xlr_train_with_noise.to_csv(os.path.join(out_dir,"XLR_train_noisy.txt.gz"), compression = "gzip", sep = "\t", index = False)
            ds_xlr_test_with_noise.to_csv(os.path.join(out_dir,"XLR_test_noisy.txt.gz"), compression = "gzip", sep = "\t", index = False)
        
        if not skip_training:
            #nsd and sd if not skipping training
            ds_nsd = pd.concat([pd.DataFrame(output_dict["x_nsd_or"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_nsd"], name='Class')], axis=1)
            ds_sd = pd.concat([pd.DataFrame(output_dict["x_sd_or"], columns =['GC_content', 'Length', 'NRC_poolNorm'] ),
                                        pd.Series(output_dict["y_sd"], name='Class')], axis=1)
            
            ds_nsd.to_csv(os.path.join(out_dir,"NSD.txt.gz"), compression = "gzip", sep = "\t", index = False)
            ds_sd.to_csv(os.path.join(out_dir,"SD.txt.gz"), compression = "gzip", sep = "\t", index = False)
            
    return output_dict
   
    