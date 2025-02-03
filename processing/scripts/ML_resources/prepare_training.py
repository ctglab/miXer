#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 13:26:12 2022

@author: elia
"""

#function to split samples from labels
#adds gaussian noise with mean = mu and variance = sigma
import os
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

#controlling randomness
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value) #fixed ambient variable value
random.seed(seed_value) #fixed random module seed
np.random.seed(seed_value) #fixed numpy's random module seed

def add_noise_to_data(mu, sigma, clean_data):
    """
    Adds Gaussian noise to the given clean data.

    Parameters:
    mu (float): The mean of the Gaussian noise distribution.
    sigma (float): The standard deviation of the Gaussian noise distribution.
    clean_data (numpy.ndarray): The clean data to which noise will be added.

    Returns:
    numpy.ndarray: The noisy data obtained by adding Gaussian noise to the clean data.
    """
    
    noise = np.random.normal(mu, sigma, clean_data.shape)
    
    noisy_data = clean_data + noise
    
    return noisy_data

def prepare_training(data, training_columns = None, split_to_train = False, split_fraction = 0.2, add_noise = False, 
                     train_samples = None, mu = 0, sigma = 0.05, seed_value = None):
    
    """
    Prepares the training data for model training and testing on chrX.

    Parameters:
    data (pandas.DataFrame): The input dataset.
    training_columns (list, optional): The list of column names to use as training features. If None, all columns except "Class" will be used. Default is None.
    split_to_train (bool, optional): Whether to split the data into training and testing sets. Default is False.
    split_fraction (float, optional): The fraction of data to be used for testing when splitting the data. Default is 0.2.
    add_noise (bool, optional): Whether to add noise to the data. Default is False.
    train_samples (int, optional): The number of samples to be used for training when splitting the data. If None, the fraction specified by split_fraction will be used. Default is None.
    mu (float, optional): The mean of the Gaussian noise distribution when adding noise to the data. Default is 0.
    sigma (float, optional): The standard deviation of the Gaussian noise distribution when adding noise to the data. Default is 0.05.
    seed_value (int, optional): The random seed value for reproducible results. Default is None.

    Returns:
    dict or tuple: If split_to_train is True, a dictionary containing the prepared training data in the following format:
                  {
                      "x_train_with_noise": numpy.ndarray,
                      "x_train_noisy_only": numpy.ndarray,
                      "x_train_clean_only": numpy.ndarray,
                      "x_test_with_noise": numpy.ndarray,
                      "x_test_noisy_only": numpy.ndarray,
                      "x_test_clean_only": numpy.ndarray,
                      "y_train_with_noise": numpy.ndarray,
                      "y_train_noisy_only": numpy.ndarray,
                      "y_train_clean_only": numpy.ndarray,
                      "y_test_with_noise": numpy.ndarray,
                      "y_test_noisy_only": numpy.ndarray,
                      "y_test_clean_only": numpy.ndarray
                  }
                  If split_to_train is False, a tuple containing the training features numpy.ndarray and the target labels numpy.ndarray.
    """
    
    dataset = data.copy(deep = True)
    if training_columns is None:
        orig_train_columns = [x for x in dataset.columns if x != "Class"]
    else:
        orig_train_columns = training_columns
        
    if split_to_train:
            
        if add_noise:
            
            labels= dataset.loc[:, dataset.columns == "Class"]
            samples = dataset.loc[:, dataset.columns != "Class"]
            
            noisy_samples = add_noise_to_data(mu, sigma, samples)
            noisy_labels = labels.copy()
            
            samples["noisy"] = False
            noisy_samples["noisy"] = True
            
            samples = pd.concat([samples, noisy_samples])
            labels = pd.concat([labels, noisy_labels])
            
            samples.reset_index(drop = True,inplace = True)
            labels.reset_index(drop=True, inplace = True)
            
            labels = labels.astype("int")
    
            #split training-test
            if train_samples == None:
                x_train, x_test, y_train, y_test = train_test_split(samples,
                                                                    labels, 
                                                                    test_size = split_fraction,
                                                                    stratify = labels,
                                                                    random_state = seed_value)
            else:
                x_train, x_test, y_train, y_test = train_test_split(samples,
                                                                    labels, 
                                                                    train_size = train_samples,
                                                                    stratify = labels,
                                                                    random_state = seed_value)
            
            #creating noisy version of x_train
            x_train_noisy_only = x_train.loc[x_train["noisy"] == True]
            x_train_clean_only = x_train.loc[x_train["noisy"] == False]
            x_test_noisy_only = x_test.loc[x_test["noisy"] == True]
            x_test_clean_only = x_test.loc[x_test["noisy"] == False]
            
            y_train_noisy_only = y_train[y_train.index.isin(x_train_noisy_only.index)]
            y_train_clean_only = y_train[y_train.index.isin(x_train_clean_only.index)]
            y_test_noisy_only = y_test[y_test.index.isin(x_test_noisy_only.index)]
            y_test_clean_only = y_test[y_test.index.isin(x_test_clean_only.index)]
            
            x_train = x_train[orig_train_columns]
            x_test = x_test[orig_train_columns]
            x_train_noisy_only = x_train_noisy_only[orig_train_columns]
            x_train_clean_only = x_train_clean_only[orig_train_columns]
            x_test_noisy_only = x_test_noisy_only[orig_train_columns]
            x_test_clean_only = x_test_clean_only[orig_train_columns]
            
            #output a dictionary, too many outputs
            output_dict = {
                "x_train_with_noise": x_train.to_numpy(),
                "x_train_noisy_only": x_train_noisy_only.to_numpy(),
                "x_train_clean_only": x_train_clean_only.to_numpy(),
                "x_test_with_noise": x_test.to_numpy(),
                "x_test_noisy_only": x_test_noisy_only.to_numpy(),
                "x_test_clean_only": x_test_clean_only.to_numpy(),
                "y_train_with_noise": np.ravel(y_train.to_numpy()),
                "y_train_noisy_only": np.ravel(y_train_noisy_only.to_numpy()),
                "y_train_clean_only": np.ravel(y_train_clean_only.to_numpy()),
                "y_test_with_noise": np.ravel(y_test.to_numpy()),
                "y_test_noisy_only": np.ravel(y_test_noisy_only.to_numpy()),
                "y_test_clean_only": np.ravel(y_test_clean_only.to_numpy())
                }
        else:
            #split training-test
            if train_samples == None:
                x_train, x_test, y_train, y_test = train_test_split(dataset[[x for x in dataset.columns if x != "Class"]],
                                                                    dataset["Class"], 
                                                                    test_size = split_fraction,
                                                                    stratify = dataset["Class"],
                                                                    random_state = seed_value)
            else:
                x_train, x_test, y_train, y_test = train_test_split(dataset[[x for x in  dataset.columns if x != "Class"]],
                                                                    dataset["Class"], 
                                                                    train_size = train_samples,
                                                                    stratify = dataset["Class"],
                                                                    random_state = seed_value)
            output_dict = {
                "x_train_clean_only": x_train.to_numpy(),
                "x_test_clean_only": x_test.to_numpy(),
                "y_train_clean_only": y_train.to_numpy(),
                "y_test_clean_only": y_test.to_numpy()
                }
            
        return output_dict
    
    else:
        return dataset[orig_train_columns].to_numpy(), dataset["Class"].to_numpy()