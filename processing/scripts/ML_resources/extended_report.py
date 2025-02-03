#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 22:29:35 2022

@author: elia
"""

#function to print a report containing useful informations about the model's performance on training datasets
import os
#import numpy as np

from ML_resources.get_min1_0_1 import get_min1_0_1
from ML_resources.misclassified_gld_nsd_sd import misclassified_gld_nsd_sd
from ML_resources.fp_fn_gld_nsd_sd import fp_fn_gld_nsd_sd
#from training_model_scores import training_model_scores

def extended_report(modelname, 
                    y_train_golden, y_test_golden, y_nsd, y_sd,
                    res_gld_train, res_gld_test, res_nsd, res_sd, savepath,
                    sc_metrics_all,
                    avg = "macro"):
    """
    Generate an extended performance report for a classification model.

    The function prints and saves an extended performance report for a classification model,
    including information about the dataset samples, label distribution, misclassified samples,
    false positives, false negatives, and evaluation scores.

    Args:
        modelname (str): The name of the model.
        y_train_golden (array-like): The true labels for the training set.
        y_test_golden (array-like): The true labels for the test set.
        y_nsd (array-like): The true labels for the NSD set.
        y_sd (array-like): The true labels for the SD set.
        res_gld_train (tuple): A tuple containing the prediction results for the training set.
            The tuple should have the following structure:
            (predictions, probabilities, accuracy, precision, recall, f1_score)
        res_gld_test (tuple): A tuple containing the prediction results for the test set.
            The tuple should have the same structure as `res_gld_train`.
        res_nsd (tuple): A tuple containing the prediction results for the NSD set.
            The tuple should have the same structure as `res_gld_train`.
        res_sd (tuple): A tuple containing the prediction results for the SD set.
            The tuple should have the same structure as `res_gld_train`.
        savepath (str): The path to save the report file.
        avg (str, optional): The averaging strategy to use for precision, recall, and F1-score.
            Supported options are 'macro' (default), 'micro', 'weighted', and None.

    Returns:
        None
    """
    
    #quanti sample nei dataset:
    gold_train_samples = len(y_train_golden)
    gold_test_samples = len(y_test_golden)
    nsd_samples = len(y_nsd)
    sd_samples = len(y_sd)
    
    #numero di -2, -1, 0 e 1
    mun_min2_gld_train, mun_min1_gld_train, num_0_gld_train, num_1_gld_train, num_2_gld_train= get_min1_0_1(y_train_golden)
    num_min2_gld_test, num_min1_gld_test, num_0_gld_test, num_1_gld_test, num_2_gld_test = get_min1_0_1(y_test_golden)
    num_min2_nsd, num_min1_nsd, num_0_nsd, num_1_nsd, num_2_nsd = get_min1_0_1(y_nsd)
    num_min2_sd, num_min1_sd, num_0_sd, num_1_sd, num_2_sd  =get_min1_0_1(y_sd)
    
    #numero di misclassificati per dataset, threshold standard
    
    mis_gld_train, mis_gld_test, mis_nsd, mis_sd  = misclassified_gld_nsd_sd(res_gld_train[0], res_gld_test[0],
                                                                              res_nsd[0], res_sd[0],
                                                                              y_train_golden, y_test_golden,
                                                                              y_nsd, y_sd)
    
    
    #numero di falsi positivi e falsi negativi, threshold standard
    result = fp_fn_gld_nsd_sd(res_gld_train[0], res_gld_test[0],
                          res_nsd[0], res_sd[0],
                          y_train_golden, y_test_golden,
                          y_nsd, y_sd)

    #spacchetta i risultati
    fp_gold_train, fn_gold_train, tpr_gold_train, fpr_gold_train, tnr_gold_train, fnr_gold_train = result[0]
    fp_gold_test, fn_gold_test, tpr_gold_test, fpr_gold_test, tnr_gold_test, fnr_gold_test = result[1]
    fp_nsd, fn_nsd, tpr_nsd, fpr_nsd, tnr_nsd, fnr_nsd = result[2]
    fp_sd, fn_sd, tpr_sd, fpr_sd, tnr_sd, fnr_sd = result[3]
    
    
    print("\n{} model performance report".format(modelname))
    print("XLR train samples: {} | XLR test samples: {} | NSD samples: {} | SD samples: {}".format(gold_train_samples, gold_test_samples, nsd_samples, sd_samples))
    print("XLR train sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}".format(mun_min2_gld_train,mun_min1_gld_train, num_0_gld_train, num_1_gld_train, num_2_gld_train))
    print("XLR test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}".format(num_min2_gld_test, num_min1_gld_test, num_0_gld_test, num_1_gld_test, num_2_gld_test))
    print("NSD test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}".format(num_min2_nsd, num_min1_nsd, num_0_nsd, num_1_nsd, num_2_nsd))
    print("SD test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}\n".format(num_min2_sd, num_min1_sd, num_0_sd, num_1_sd, num_2_sd))
    #print("\nThreshold standard ( x = 1 (or -1) if prediction confidence >= 0.5)")
    print("Total misclassified XLR train samples: {}".format(mis_gld_train))
    print("Total misclassified XLR test samples: {}".format(mis_gld_test))
    print("Total misclassified NSD samples: {}".format(mis_nsd))
    print("Total misclassified SD samples: {}".format(mis_sd))
    print("XLR train false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%".format(fp_gold_train, fn_gold_train,tpr_gold_train*100, fpr_gold_train*100, tnr_gold_train*100, fnr_gold_train*100 ))
    print("XLR test false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%".format(fp_gold_test, fn_gold_test,tpr_gold_test*100, fpr_gold_test*100, tnr_gold_test*100, fnr_gold_test*100 ))
    print("NSD false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%".format(fp_nsd, fn_nsd,tpr_nsd*100, fpr_nsd*100, tnr_nsd*100, fnr_nsd*100))
    print("SD false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%".format(fp_sd, fn_sd, tpr_sd*100, fpr_sd*100, tnr_sd*100, fnr_sd*100))
    print("XLR train Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}".format(res_gld_train[2]*100,res_gld_train[3]*100,res_gld_train[4]*100,res_gld_train[5]*100))
    print("XLR test Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}".format(res_gld_test[2]*100,res_gld_test[3]*100,res_gld_test[4]*100,res_gld_test[5]*100))
    print("NSD Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}".format(res_nsd[2]*100,res_nsd[3]*100,res_nsd[4]*100,res_nsd[5]*100))
    print("SD Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}\n".format(res_sd[2]*100,res_sd[3]*100,res_sd[4]*100,res_sd[5]*100))
    
    print("################# Single class metrics #################")
    for key in sc_metrics_all.keys():
        print(f"Single class metrics for dataset: {key}")
        # Print row index and all column values
        metrics_df = sc_metrics_all[key].multiply(100)
        for idx, row in metrics_df.iterrows():
            print(f"Class: {idx}")
            row_str = ', '.join([f"{key}: {value:.2f}%" for key, value in row.items()])
            print(row_str)
            print("\n\n")
    
    #salva su file
    txtname = "{}_model_perf_report.txt".format(modelname)
    f = open(os.path.join(savepath, txtname), "w+")
    f.write("XLR train samples: {} | XLR test samples: {} | NSD samples: {} | SD samples: {}\n".format(gold_train_samples, gold_test_samples, nsd_samples, sd_samples))
    f.write("XLR train sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}\n".format(mun_min2_gld_train,mun_min1_gld_train, num_0_gld_train, num_1_gld_train, num_2_gld_train))
    f.write("XLR test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}\n".format(num_min2_gld_test, num_min1_gld_test, num_0_gld_test, num_1_gld_test, num_2_gld_test))
    f.write("NSD test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}\n".format(num_min2_nsd, num_min1_nsd, num_0_nsd, num_1_nsd, num_2_nsd))
    f.write("SD test sample labels | -2: {} | -1: {} | 0: {} | 1: {} | 2+:{}\n".format(num_min2_sd, num_min1_sd, num_0_sd, num_1_sd, num_2_sd))
    #f.write("\nThreshold standard ( x = 1 (or -1) if prediction confidence >= 0.5)\n")
    f.write("Total misclassified XLR train samples: {}\n".format(mis_gld_train))
    f.write("Total misclassified XLR test samples: {}\n".format(mis_gld_test))
    f.write("Total misclassified NSD samples: {}\n".format(mis_nsd))
    f.write("Total misclassified SD samples: {}\n".format(mis_sd))
    f.write("XLR train false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%\n".format(fp_gold_train, fn_gold_train,tpr_gold_train*100, fpr_gold_train*100, tnr_gold_train*100, fnr_gold_train*100 ))
    f.write("XLR test false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%\n".format(fp_gold_test, fn_gold_test,tpr_gold_test*100, fpr_gold_test*100, tnr_gold_test*100, fnr_gold_test*100 ))
    f.write("NSD false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%\n".format(fp_nsd, fn_nsd,tpr_nsd*100, fpr_nsd*100, tnr_nsd*100, fnr_nsd*100 ))
    f.write("SD false positives: {} | false negatives: {} | TPR: {:.2f}% | FPR: {:.2f}% | TNR: {:.2f}% | FNR: {:.2f}%\n".format(fp_sd, fn_sd, tpr_sd*100, fpr_sd*100, tnr_sd*100, fnr_sd*100))
    f.write("XLR train Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}\n".format(res_gld_train[2]*100,res_gld_train[3]*100,res_gld_train[4]*100,res_gld_train[5]*100))
    f.write("XLR test Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}\n".format(res_gld_test[2]*100,res_gld_test[3]*100,res_gld_test[4]*100,res_gld_test[5]*100))
    f.write("NSD Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}\n".format(res_nsd[2]*100,res_nsd[3]*100,res_nsd[4]*100,res_nsd[5]*100))
    f.write("SD Accuracy: {:.2f} | Precision: {:.2f} | Recall: {:.2f} | F1-score: {:.2f}\n".format(res_sd[2]*100,res_sd[3]*100,res_sd[4]*100,res_sd[5]*100))
    f.write("################# Single class metrics #################\n")
    for key in sc_metrics_all.keys():
        f.write(f"Single class metrics for dataset: {key}\n")
        metrics_df = sc_metrics_all[key].multiply(100)
        for idx, row in metrics_df.iterrows():
            f.write(f"Class: {idx}\n")
            row_str = ', '.join([f"{key}: {value:.2f}%" for key, value in row.items()])
            f.write(row_str)
            f.write("\n")
    f.close()
    
    return
