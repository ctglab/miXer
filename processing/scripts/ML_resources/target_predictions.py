#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 09:31:34 2022

@author: elia
"""

#function to make predictions on test samples (like NA12878)
import os
import pandas as pd
import numpy as np


def target_predictions(target, target_name, modelname, training_columns, 
                       scaler, clf, output_folder, dataset_dir_name,
                       orig_columns,
                       force_mnorm = True, verbose = 1, skip_if_present = False):
    """
    Performs target predictions using a machine learning model.

    Parameters:
    target (pandas.DataFrame): The target data to be predicted.
    target_name (str): The name of the target.
    modelname (str): The name of the machine learning model.
    training_columns (list): The list of columns used for training the model.
    scaler (object): The scaler object used for scaling the data.
    clf (object): The machine learning classifier model.
    output_folder (str): The path to the output folder where predictions will be saved.
    dataset_dir_name (str): The name of the dataset directory.
    orig_columns: The original columns of the target data.
    force_mnorm (bool, optional): Whether to force median normalization. Defaults to True.
    verbose (int, optional): The verbosity level. Defaults to 1.
    skip_if_present (bool, optional): Whether to skip if predictions are already present. Defaults to False.

    Returns:
    None
    """
    
    target_copy = target.copy()
    if "TARGET" not in target_name:
        basename = target_name + "_TARGET_{}_pred".format(modelname)
    else:
        basename = target_name + "_{}_pred".format(modelname)
            
    savepath = os.path.join(output_folder, target_name)
    target_out = os.path.join(savepath, modelname + "_" + dataset_dir_name)
    
    if os.path.exists(savepath) == False:
        os.makedirs(savepath)
    
    if os.path.exists(target_out) == False:
        os.makedirs(target_out)
        
    exec_checks = [False]
    
    savename = basename + ".txt.gz"
    if skip_if_present:
        exec_checks.append(os.path.isfile(os.path.join(target_out, savename)))
    
    f = open(os.path.join(target_out, "{}_MLcall_log.txt".format(basename)), "w+")
    f.write("Model ID: {}\n".format(modelname))
    
    if not all(exec_checks):
        
        #forcing median normalization of NRC poolnorm, always done, no more a parameter
        nrc_median = target_copy[~target_copy["Chr"].isin(["chrX", "ChrX", "X"])]["NRC_poolNorm"].quantile(0.5)
        
        if verbose >= 1:
            print("Forcing median normalization for sample {}. Old median {:.2f}".format(target_name,nrc_median))
            f.write("Forcing median normalization for sample {}. Old median {:.2f}".format(target_name,nrc_median))
        
        target_copy['NRC_poolNorm'] = target_copy['NRC_poolNorm'] - nrc_median
        if verbose > 2:
            print("New autosome median: {}".format(target_copy[~target_copy["Chr"].isin(["chrX", "ChrX", "X"])]["NRC_poolNorm"]. quantile(0.5)))
            print("ChrX median: {}".format(target_copy[target_copy["Chr"].isin(["chrX", "ChrX", "X"])]["NRC_poolNorm"]. quantile(0.5)))
    
        x_target = target_copy[training_columns].values
            
        if scaler != None:
            x_target_scaled = scaler.transform(x_target)
        else:
            x_target_scaled = x_target
        
        #making predictions with ML-model
        if verbose > 2:
            print("Making predictions")
        y_pred = clf.predict(x_target_scaled)
        y_pred_proba = clf.predict_proba(x_target_scaled)
        y_pred_prediction_probs = np.amax(y_pred_proba, axis=1)
            
        predictions_column_name = "{}_pred".format(modelname)
        
        base_pred_proba_name = "{}_pred_proba".format(modelname)
        confidence_column_name = "{}_pred_confidence".format(modelname)
        
        pred_proba_min_two = base_pred_proba_name + "_(-2)"
        pred_proba_min_one = base_pred_proba_name + "_(-1)"
        pred_proba_zero = base_pred_proba_name + "_(0)"
        pred_proba_one = base_pred_proba_name + "_(1)"
        pred_proba_two = base_pred_proba_name + "_(2)"
        
        #adding predictions and prediction probabilities to target data
        target_copy = target_copy.join(pd.DataFrame(y_pred, columns = [predictions_column_name]))
        target_copy = target_copy.join(pd.DataFrame(y_pred_prediction_probs, columns = [base_pred_proba_name]))
        
        #calculate mappability mediated SVM confidence score
        #enforce numerical dtype to avoid typerror
        target_copy[base_pred_proba_name] = pd.to_numeric(target_copy[base_pred_proba_name], errors='coerce')
        target_copy[confidence_column_name] = target_copy[base_pred_proba_name] * target["Mappability"]
        
        #append prediction probability columns to target dataframe
        if y_pred_proba.shape[1] == 3:
            target_copy = target_copy.join(pd.DataFrame(y_pred_proba, columns = [pred_proba_min_one, pred_proba_zero, pred_proba_one]))
        elif y_pred_proba.shape[1] == 4:
            target_copy = target_copy.join(pd.DataFrame(y_pred_proba, columns = [pred_proba_min_two, pred_proba_min_one, pred_proba_zero, pred_proba_one]))
        elif y_pred_proba.shape[1] == 5:
            target_copy = target_copy.join(pd.DataFrame(y_pred_proba, columns = [pred_proba_min_two, pred_proba_min_one, pred_proba_zero, pred_proba_one, pred_proba_two]))
             
        #saving predictions and filtered predictions
        if verbose >0:
            print("Calls for sample {}: {}".format(target_name, target_copy.loc[target_copy[predictions_column_name] != 0].shape[0]))
        target_copy.to_csv(os.path.join(target_out,savename), compression = "gzip", sep = "\t", index = False)
        
        ddel_calls = target_copy.loc[target_copy[predictions_column_name] == -2].shape[0]
        del_calls = target_copy.loc[target_copy[predictions_column_name] == -1].shape[0]
        wt_calls = target_copy.loc[target_copy[predictions_column_name] == 0].shape[0]
        dup_calls = target_copy.loc[target_copy[predictions_column_name] == 1].shape[0]
        mdup_calls = target_copy.loc[target_copy[predictions_column_name] == 2].shape[0]
        
        f.write("Model calls: DDELs: {} | DELs: {} | WTs: {} | DUPs: {} | MDUPs: {}\n".format(ddel_calls, del_calls, wt_calls, dup_calls, mdup_calls))
        
        f.write("Predicted as DDEL median NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -2]["NRC_poolNorm"].quantile(0.5)))
        f.write("Predicted as DDEL min NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -2]["NRC_poolNorm"].min()))
        f.write("Predicted as DDEL max NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -2]["NRC_poolNorm"].max()))
        f.write("Predicted as DDEL 1st Q NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -2]["NRC_poolNorm"].quantile(0.25) ))
        f.write("Predicted as DDEL 3rd Q NRC: {}\n\n".format(target_copy.loc[target_copy[predictions_column_name] == -2]["NRC_poolNorm"].quantile(0.75) ))
        
        f.write("Predicted as DEL median NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -1]["NRC_poolNorm"].quantile(0.5)))
        f.write("Predicted as DEL min NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -1]["NRC_poolNorm"].min()))
        f.write("Predicted as DEL max NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -1]["NRC_poolNorm"].max()))
        f.write("Predicted as DEL 1st Q NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == -1]["NRC_poolNorm"].quantile(0.25) ))
        f.write("Predicted as DEL 3rd Q NRC: {}\n\n".format(target_copy.loc[target_copy[predictions_column_name] == -1]["NRC_poolNorm"].quantile(0.75) ))
        
        
        f.write("Predicted as DUP median NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 1]["NRC_poolNorm"].quantile(0.5)))
        f.write("Predicted as DUP min NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 1]["NRC_poolNorm"].min()))
        f.write("Predicted as DUP max NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 1]["NRC_poolNorm"].max()))
        f.write("Predicted as DUP 1st Q NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 1]["NRC_poolNorm"].quantile(0.25) ))
        f.write("Predicted as DUP 3rd Q NRC: {}\n\n".format(target_copy.loc[target_copy[predictions_column_name] == 1]["NRC_poolNorm"].quantile(0.75) ))
        
        f.write("Predicted as MDUP median NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 2]["NRC_poolNorm"].quantile(0.5)))
        f.write("Predicted as MDUP min NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 2]["NRC_poolNorm"].min()))
        f.write("Predicted as MDUP max NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 2]["NRC_poolNorm"].max()))
        f.write("Predicted as MDUP 1st Q NRC: {}\n".format(target_copy.loc[target_copy[predictions_column_name] == 2]["NRC_poolNorm"].quantile(0.25) ))
        f.write("Predicted as MDUP 3rd Q NRC: {}\n\n".format(target_copy.loc[target_copy[predictions_column_name] == 2]["NRC_poolNorm"].quantile(0.75) ))
        
        f.close()
        
        del target_copy
        del x_target
        del x_target_scaled
        
        del y_pred
        del y_pred_proba
    else:
        if verbose > 2:
            print("Sample {} already present, skipping".format(target_name))
        
    return