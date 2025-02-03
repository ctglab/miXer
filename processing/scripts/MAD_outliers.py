#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 09:03:34 2023

@author: elia
"""

#MAD calculation script for outlier detection

import os
import argparse

import pandas as pd

#Defining default values
DEFAULT_MAD_THR = 0.1

#parsing arguments
ap = argparse.ArgumentParser(description = "Test for MADs!")

ap.add_argument('-exca', "--excavator_windows", type=str, help='Excavator2 windows folder path')
ap.add_argument("-thr", "--mad_threshold", default = DEFAULT_MAD_THR, type = float, help= "SegMean MAD threshold to consider a sample an outlier. Default = {}".format(DEFAULT_MAD_THR))
ap.add_argument("-ename", "--expname", type = str, default = None, help = "Experiment name (will be used as output folder name if not None)")
ap.add_argument("-outpath", "--output_folder", type = str, default = os.getcwd(), help = "Main output dir folder. Default = {}".format(os.getcwd()))

args = vars(ap.parse_args())
exca_folder = args["excavator_windows"]
mad_thr = args["mad_threshold"]
expname = args["expname"]
outpath = args["output_folder"]

#creating output folder
if expname is not None:
    output_directory = os.path.join(outpath, expname)
else:
    output_directory = os.path.join(outpath, "Excavator2_MAD_Outlier_Detection")

if not os.path.isdir(output_directory):
    os.mkdir(output_directory)

#get files list in exca folder path
exca_subfolders = [f for f in os.listdir(exca_folder)]

segmean_mads = []

print("Calculating SegMean signal variation metrics")
#cycle over all the samples present in both callers
for samplename in exca_subfolders:
    
    #get path for exca fastcall file
    exca_hslm_path = os.path.join(exca_folder, os.path.join(samplename, "HSLMResults_" + samplename +".txt"))
    
    #Load excavator Fastcall
    exca_hslm = pd.read_table(exca_hslm_path)
    
    segmean_mads.append(exca_hslm["SegMean"].loc[exca_hslm["Chromosome"] != "chrX"].mad())
    
metrics = {"samplename": exca_subfolders,
           "segmean_mads": segmean_mads}

#Saving all samples and their MADs
metrics_df = pd.DataFrame(data=metrics)
metrics_df.to_csv(os.path.join(output_directory, "All_Samples_with_MAD.csv"), sep = "\t", index=False)

#Filtering and keeping only files that have MADs under the threshold
to_keep = metrics_df.loc[metrics_df["segmean_mads"] <= mad_thr]
to_keep.to_csv(os.path.join(output_directory, "Samples_ToKeep_with_MAD.csv"), sep = "\t", index=False)


