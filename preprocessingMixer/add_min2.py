#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 10:57:38 2022

@author: elia
"""

#script that adds simulated min2 data to a dataset (XLR, NSD or SD)
#adds the right amount of min2 (1/3 of the size of the dataset) to keep having a balanced dataset

import pandas as pd
import os
import argparse
import numpy as np

def check_none(string):
        if string.lower() in ["none"]:
            return None
        else:
            return string
            
ap = argparse.ArgumentParser()
ap.add_argument("-ddir", "--dataset_folder", required = True, help = "Path to dataset folder.")
ap.add_argument("-did", "--dataset_id", required = True, help = "Dataset ID. Either XLR, NSD or SD.")
ap.add_argument("-out_old", "--old_data_output_path", required = True, help = "Path to output folder for old dataset.\nThe script will remove the old dataset from the folder, this is just to save the old version somewhere. Use 'none' to overwrite only")
ap.add_argument("-min2", "--min2_dataset_path", required = True, help = "Path to min2 dataset.")

args = vars(ap.parse_args())
data_fold = args["dataset_folder"]
data_id = args["dataset_id"]
min2_path = args["min2_dataset_path"]

out_old = check_none(args["old_data_output_path"])

allowed_ids = ["xlr", "nsd", "sd", "nosegdup", "segdup"]

if data_id.lower() in allowed_ids:
    check_id = True
else:
    check_id = False
    
#format dataset id for filename matching
if data_id.lower() in ["xlr", "xl", "xr", "x"]:
    data_id = "XLR"
elif data_id.lower() in ["nsd", "nosegdup", "noseg"]:
    data_id = "noSegDup"
elif data_id.lower() in ["sd", "segdup", "seg"]:
    data_id = "SegDup"

seed_value = 42

#define columns 
columns = ["GC_content", "Length",   "MeanCvg",  "NRC_poolNorm", "NRC_delta", "Class"]

#check if dataset id is correct
if check_id:
    
    #TODO from data id get filename
    #read filenames in data folder
    f = []
    for (dirpath, dirnames, filenames) in os.walk(data_fold):
        f.extend(filenames)
        break
    
    #take only the filename associated to the tag
    old_data_name = [x for x in f if data_id.lower() in x.lower()][0]
    
    #split on "." to get prefix
    old_data_prefix = old_data_name.split(".")[0]
    
    if "Min2" in old_data_prefix or "DDup" in old_data_prefix:
        print("Warning! Probably you already ran the script for the {} dataset.".format(data_id))
    
    new_data_name = old_data_prefix + "_Min2.txt.gz"
    
    data_path = os.path.join(data_fold, old_data_name)
    
    #load datasets
    dataset = pd.read_table(data_path, sep = "\t")
    #if user wants, save old version
    if out_old is not None:
        if not os.path.isdir(out_old):
            os.mkdir(out_old)
            
        dataset.to_csv(os.path.join(out_old, old_data_name), sep = "\t", compression ="gzip", index = False)
    
    #remove old dataset
    os.remove(data_path)
    
    min2_dataset = pd.read_table(min2_path, sep = "\t")
    
    #check if min2 data already present, if so remove it and substitute it with new min2 data 
    if -2 in dataset["Class"].unique():
        
        dataset = dataset.loc[dataset["Class"] != -2]
    
    #check the number of deletions in dataset (it should be equal to the number of duplications/wts)
    num_del = dataset.loc[dataset["Class"] == -1].shape[0]
    
    #subsample to keep balance between classes
    if num_del < min2_dataset.shape[0]:
        #There are less deletions than min2 data available
        #subsampling min2 dataset to keep class balance in dataset
        min2_dataset = min2_dataset.sample(n = num_del, ignore_index = True, random_state = seed_value)
    
    dataset = pd.concat([dataset, min2_dataset])
    
    #shuffling
    dataset = dataset.sample(frac=1, random_state = seed_value).reset_index(drop=True)
    
    #removing useless columns
    dataset = dataset[columns]
    
    #saving in output folder, in place of the old XLR
    dataset.to_csv(os.path.join(data_fold, new_data_name), sep = "\t", compression ="gzip", index = False)
    
else:
    print("Wrong dataset ID supplied. Execution aborted.")









