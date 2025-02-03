#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 15:37:40 2021
​
@author: elia
"""

#merged dataset maker

#TODO descrizione

import pandas as pd
import numpy as np
import os
import random
import argparse
import sys
import time
import re
import datetime

timestamp = time.gmtime()

ap = argparse.ArgumentParser()
ap.add_argument("-ddir", "--datasets_directory", required = True, help = "Path to directory containing each kit's XLR, NoSegDup and SegDup dataset.")

args = vars(ap.parse_args())
data_dir = args["datasets_directory"]
curr_dir = os.getcwd()

dirname = os.path.basename(os.path.dirname(data_dir))
merged_dir = os.path.join(data_dir, "merged")

########## Generating merged dataset

#controlling randomness
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value) #fixed ambient variable value
random.seed(seed_value) #fixed random module seed
np.random.seed(seed_value) #fixed numpy's random module seed

#logging the execution
logfile_name = "merged_maker_dataset_{}_logfile_{}.txt".format(dirname, time.strftime("%c", timestamp))
log = open(os.path.join(curr_dir, logfile_name), "w+")

log.write("Log for merged_maker.py execution\nMain folder directory: {}\n".format(data_dir))
log.write("Executed on {} Timezone: {}\n".format(time.strftime("%c", timestamp), time.strftime("%Z", timestamp))) #%c for locale’s appropriate date and time representation using time.strftime function

#check if merged folder is present, to save it
if os.path.exists(merged_dir) == False:
    #check if dataset folder(s) are present
    print("Merged datasets folder not found. Creating...")
    log.write("Merged datasets folder not found. Creating...\n")
    #TODO delete old if present?
    os.makedirs(merged_dir)
    
if os.path.exists(merged_dir) == False:
    os.makedirs(merged_dir)

#check if dataset folder is not empty (script already ran)

if len(os.listdir(merged_dir)) != 0:
    print("Non-empty Merged datasets directory detected")
    log.write("Non-empty Merged datasets directory detected\n")
    
    inp = str(input("Do you want to continue the execution? Old datasets will be overwritten! "))
    log.write("Do you want to continue the execution? Old datasets will be overwritten! " + inp + "\n")
    
    if inp.lower() in ["y", "yes", "yeah", "ok", "1", "yup", "true", "t"]:
        print("Moving on with execution")
        log.write("Moving on with execution\n")
    elif inp.lower() in ["no", "n", "nope", "0", "false", "f"]:
        log.write("Exiting...")
        log.close()
        sys.exit("Exiting...")
    else:
        log.write("Invalid input. Exiting...")
        log.close()
        sys.exit("Invalid input. Exiting...")
            
#cycle over the dataset folders, excluding merged

#takes the list of the folders in data directory --> the datasets
dataset_list = [] 
for (dirpath, dirnames, filenames) in os.walk(data_dir):
    dataset_list.extend(dirnames)
    break


#remove merged if present, must not consider it, old will be overwritten
try:
    dataset_list.remove("merged")
except:
    pass

#must check which dataset is the minimal one
min_exons_gld = 1e10
min_exons_nsd = 1e10
min_exons_sd = 1e10

print("Checking sizes of datasets")
log.write("Folder list: {}\n".format(dataset_list))
log.write("Checking sizes of datasets\n")

#TODO log info on datasets
for kit in dataset_list:
    
    #####TODO controllare, a me non funziona (Elia)
    # kit_folder = os.path.join(data_dir, kit)
    # subfolders = [ os.path.basename(f.path) for f in os.scandir(kit_folder) if f.is_dir() ]
    # matchDir = re.compile('.*datasets_[0-9-]')
    # filtered = [folder for folder in subfolders if matchDir.match(folder)]
    # if not filtered:
    #    continue
    # else:
    #    dataset_folder = sorted(filtered, reverse=True)[0]       
    #    full_path= os.path.join(kit_folder,dataset_folder)
    
    #hotfix
    full_path = os.path.join(data_dir, kit)
    
    #load the datasets
    #check if they exist
    gld_path = os.path.join(full_path, "ALL_SAMPLE_XLR.txt.gz")
    nsd_path = os.path.join(full_path, "ALL_SAMPLE_noSegDup.txt.gz")
    sd_path = os.path.join(full_path, "ALL_SAMPLE_SegDup.txt.gz")
    
    gld_ok = os.path.isfile(gld_path)
    nsd_ok = os.path.isfile(nsd_path)
    sd_ok = os.path.isfile(sd_path)
    
    data_ok = all([gld_ok, nsd_ok, sd_ok]) #returns true if all items in iterable are true
   
    if data_ok == True:
        logname = "{}_dataset_info.txt".format(kit)
        f = open(os.path.join(full_path, logname), "w+")
        
        gld = pd.read_csv(os.path.join(full_path, "ALL_SAMPLE_XLR.txt.gz"), compression = "gzip", sep = "\t")
        nsd = pd.read_csv(os.path.join(full_path, "ALL_SAMPLE_noSegDup.txt.gz"), compression = "gzip", sep = "\t")
        sd = pd.read_csv(os.path.join(full_path, "ALL_SAMPLE_SegDup.txt.gz"), compression = "gzip", sep = "\t")
        
        gld_ex = gld.shape[0]
        nsd_ex = nsd.shape[0]
        sd_ex = sd.shape[0]
        
        gld_dels = gld.loc[gld["Class"]== -1].shape[0]
        gld_wt = gld.loc[gld["Class"]== 0].shape[0]
        gld_dups = gld.loc[gld["Class"]== 1].shape[0]
        
        nsd_dels = nsd.loc[nsd["Class"] == -1].shape[0]
        nsd_wt = nsd.loc[nsd["Class"] == 0].shape[0]
        nsd_dups = nsd.loc[nsd["Class"] == 1].shape[0]
        
        sd_dels = sd.loc[sd["Class"] == -1].shape[0]
        sd_wt = sd.loc[sd["Class"] == 0].shape[0]
        sd_dups = sd.loc[sd["Class"] == 1].shape[0]
        
        f.write("XLR dataset exons: {} | Dels: {} | WT: {} | Dups: {}\n".format(gld_ex,gld_dels, gld_wt, gld_dups))
        f.write("NoSegdup dataset exons: {} | Dels: {} | WT: {} | Dups: {}\n".format(nsd_ex,nsd_dels, nsd_wt, nsd_dups))
        f.write("SegDup dataset exons: {} | Dels: {} | WT: {} | Dups: {}\n".format(sd_ex,sd_dels, sd_wt, sd_dups))
        
        f.close()
        
        if gld_ex < min_exons_gld:
            min_exons_gld = gld_ex
            min_gld = kit
        if nsd_ex < min_exons_nsd:
            min_exons_nsd = nsd_ex
            min_nsd = kit
        if sd_ex < min_exons_sd:
            min_exons_sd = sd_ex
            min_sd = kit
    else:
        print("Dataset from {} kit not found, moving on".format(kit))
        log.write("Dataset from {} kit not found, moving on\n".format(kit))
        
print("Merged dataset will be created with {} samples: XLR: {} | noSegDup: {} | SegDup: {}\nThe number of exons matches that of the dataset {} which is the smallest.\n".format(min_exons_gld+min_exons_nsd+min_exons_sd, min_exons_gld, min_exons_nsd, min_exons_sd,min_gld))
log.write("Merged dataset will be created with {} samples: XLR: {} | noSegDup: {} | SegDup: {}\nThe number of exons matches that of the dataset {} which is the smallest.\n".format(min_exons_gld+min_exons_nsd+min_exons_sd, min_exons_gld, min_exons_nsd, min_exons_sd,min_gld))

merged_gold = []
merged_nsd = []
merged_sd = []

sampled_merged_gold_train = []
sampled_merged_nsd_train = []
sampled_merged_sd_train = []

sampled_merged_gold_test = []
sampled_merged_nsd_test = []
sampled_merged_sd_test = []

#fraction of exons to get from each dataset to build the merged datasets
gold_fraction = int(min_exons_gld/len(dataset_list))
nsd_fraction = int(min_exons_nsd/len(dataset_list))
sd_fraction = int(min_exons_sd/len(dataset_list))

#creating the merged datasets
print("Creating the merged datasets")
log.write("Creating the merged datasets\n")

for kit in dataset_list:
    
    #####TODO controllare, a me non funziona (Elia)
    # kit_folder = os.path.join(data_dir, kit)
    # subfolders = [ os.path.basename(f.path) for f in os.scandir(kit_folder) if f.is_dir() ]
    # matchDir = re.compile('.*datasets_[0-9-]')
    # filtered = [folder for folder in subfolders if matchDir.match(folder)]
    # if not filtered:
    #    continue
    # else:
    #    dataset_folder = sorted(filtered, reverse=True)[0]       
    #    full_path= os.path.join(kit_folder,dataset_folder)
    
    #hotfix
    full_path = os.path.join(data_dir, kit)
    
    #load the datasets
    #check if they exist
    gld_path = os.path.join(full_path, "ALL_SAMPLE_XLR.txt.gz")
    nsd_path = os.path.join(full_path, "ALL_SAMPLE_noSegDup.txt.gz")
    sd_path = os.path.join(full_path, "ALL_SAMPLE_SegDup.txt.gz")
    ddup_path = os.path.join(full_path, "ALL_SAMPLE_DDUP.txt.gz")
    
    gld_ok = os.path.isfile(gld_path)
    nsd_ok = os.path.isfile(nsd_path)
    sd_ok = os.path.isfile(sd_path)
    ddup_ok = os.path.isfile(ddup_path)
    
    data_ok = all([gld_ok, nsd_ok, sd_ok]) #returns true if all items in iterable are true
   
    if data_ok == True:
        gld = pd.read_csv(gld_path, compression = "gzip", sep = "\t")
        nsd = pd.read_csv(nsd_path, compression = "gzip", sep = "\t")
        sd = pd.read_csv(sd_path, compression = "gzip", sep = "\t")
                         
        gld_del = gld.loc[gld["Class"] == -1]
        gld_nc =  gld.loc[gld["Class"] == 0]
        gld_dup = gld.loc[gld["Class"] == 1]
        
        nsd_del = nsd.loc[nsd["Class"] == -1]
        nsd_nc =  nsd.loc[nsd["Class"] == 0]
        nsd_dup = nsd.loc[nsd["Class"] == 1]
        
        sd_del = sd.loc[sd["Class"] == -1]
        sd_nc =  sd.loc[sd["Class"] == 0]
        sd_dup = sd.loc[sd["Class"] == 1]
        
        #get del, nc and dup samples list
        del_samples = set(gld_del["ID"])
        nc_samples = set(gld_nc["ID"])
        dup_samples = set(gld_dup["ID"]) #IDs are the same for gld, nsd and sd datasets
        
        #adding exons to relative lists (should be faster than using df.append)
        
        #classic merged
        merged_gold.append(gld_del.sample(n = int(gold_fraction/3), random_state = seed_value).values)
        merged_gold.append(gld_nc.sample(n = int(gold_fraction/3), random_state = seed_value).values)
        merged_gold.append(gld_dup.sample(n = int(gold_fraction/3), random_state = seed_value).values)
        
        merged_nsd.append(nsd_del.sample(n = int(nsd_fraction/3), random_state = seed_value).values)
        merged_nsd.append(nsd_nc.sample(n = int(nsd_fraction/3), random_state = seed_value).values)
        merged_nsd.append(nsd_dup.sample(n = int(nsd_fraction/3), random_state = seed_value).values)
        
        merged_sd.append(sd_del.sample(n = int(sd_fraction/3), random_state = seed_value).values)
        merged_sd.append(sd_nc.sample(n = int(sd_fraction/3), random_state = seed_value).values)
        merged_sd.append(sd_dup.sample(n = int(sd_fraction/3), random_state = seed_value).values)
        
        if ddup_ok:
            log.write("Found DDUP training data for kit {}. Adding an equal number of exons to golden set.".format(kit))
            ddup = pd.read_csv(ddup_path, compression = "gzip", sep = "\t")
            if ddup.shape[0] > gld_dup.shape[0]:
                ddup = ddup.sample(n = gld_dup.shape[0], ignore_index = True, random_state = seed_value)
            merged_gold.append(ddup.values)
    else:
        print("Dataset from {} kit not found, moving on".format(kit))
        log.write("Dataset from {} kit not found, moving on\n".format(kit))
        
#creating classic merged datasets
merged_gold_df = pd.DataFrame(np.row_stack(merged_gold), columns = gld.columns)
merged_nsd_df = pd.DataFrame(np.row_stack(merged_nsd), columns = gld.columns)
merged_sd_df =pd.DataFrame(np.row_stack(merged_sd), columns = gld.columns)

#saving datasets
print("Done, saving NSD, SD,datasets")
log.write("Done, saving datasets\n")
log.close()

merged_gold_df.to_csv(os.path.join(merged_dir, "ALL_SAMPLE_XLR.txt.gz"), compression = "gzip", sep = "\t", index = False)
merged_nsd_df.to_csv(os.path.join(merged_dir, "ALL_SAMPLE_noSegDup.txt.gz"), compression = "gzip", sep = "\t",  index = False)
merged_sd_df.to_csv(os.path.join(merged_dir, "ALL_SAMPLE_SegDup.txt.gz"), compression = "gzip", sep = "\t",  index = False)