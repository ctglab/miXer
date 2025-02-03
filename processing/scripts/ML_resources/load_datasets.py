#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 12:49:05 2022

@author: elia
"""

#dataset loading function for training script

import os
import pandas as pd
import re

def load_datasets(dataset_dir):
    """
    Loads datasets from the specified directory based on file naming conventions.

    The function expects the following naming conventions for the dataset files:
    - Golden dataset (XLR exons): File name must contain 'XLR'.
    - NoSegDup dataset: File name must contain 'noSegDup'.
    - SegDup dataset: File name must contain 'SegDup'.

    Args:
        dataset_dir (str): Path to the dataset directory.

    Returns:
        tuple: Three pandas DataFrames representing the loaded datasets in the following order:
            - Golden dataset (XLR exons)
            - NoSegDup dataset
            - SegDup dataset
        If a dataset is not found or could not be loaded, the corresponding DataFrame will be None.
    """
    gld = None
    nsd = None
    sd = None

    file_patterns = {
        "XLR": re.compile(r".*XLR.*"),
        "noSegDup": re.compile(r".*noSegDup.*"),
        "SegDup": re.compile(r".*SegDup.*")
    }

    for filename in os.listdir(dataset_dir):
        filepath = os.path.join(dataset_dir, filename)

        if ".txt" in filename:
            for pattern_key, pattern in file_patterns.items():
                if pattern.match(filename):
                    if pattern_key == "XLR":
                        gld = pd.read_csv(filepath, compression="gzip", sep="\t")
                    elif pattern_key == "noSegDup":
                        nsd = pd.read_csv(filepath, compression="gzip", sep="\t")
                    elif pattern_key == "SegDup":
                        sd = pd.read_csv(filepath, compression="gzip", sep="\t")
                    break

    return gld, nsd, sd