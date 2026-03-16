#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 09:34:22 2023

@author: elia
"""
import os
import pandas as pd
from ML_resources.target_predictions import target_predictions

def parallel_predictions(folder_path, target_dfs_for_thread,
                    add_noise,
                    usecase_output_folder, model_tag, foldername, modelname,
                    training_columns, scaler, clf,
                    force_test_median_normalization,
                    train_samples, split_fraction, skip_tested):


    for target_df in target_dfs_for_thread:
        df = pd.read_csv(os.path.join(folder_path, target_df), sep="\t")

        # Rename columns for compatibility with target_predictions and HMM
        df = df.rename(columns={'chrom': 'Chr', 'NRC_poolNorm': 'NRC_poolNorm'})

        target_name = target_df.replace("_miXer_data.tsv", "")

        # Per-sample output folder (e.g., {base}/{sample_id}_SVC)
        sample_output_folder = os.path.join(usecase_output_folder, target_name + "_" + model_tag)

        if train_samples is None:
            out_foldername = foldername + "_TestFr_" + str(split_fraction)
        else:
            out_foldername = foldername + "_TrainSamp_" + str(train_samples)

        out_foldername = out_foldername + "_Noise_" + str(add_noise) + "_mnorm" + str(force_test_median_normalization)

        target_predictions(df, target_name, modelname, training_columns,
                                scaler, clf, sample_output_folder, out_foldername, training_columns,
                                force_mnorm = force_test_median_normalization,
                                skip_if_present = skip_tested)
        del df

    return
