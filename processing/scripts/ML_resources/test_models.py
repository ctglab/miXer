#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 20:21:19 2022

@author: elia
"""

#function that tests the models on XLR test set, XLR noisy test set, NSD and SD datasets
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from sklearn.metrics import confusion_matrix
import pandas as pd
import datetime
import os

from ML_resources.test_model_prediction_and_scores import test_model_prediction_and_scores
from ML_resources.extended_report import extended_report
from ML_resources.plot_confusion_matrix import plot_confusion_matrix
from ML_resources.plot_roc_curve import plot_roc_curve
from ML_resources.single_class_metrics import single_class_metrics

font = {'fontname':'Arial'}
plt.rcParams.update({'font.size': 15})
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')

def plot_performance_dotplot(modelname, savepath,
                             gld_train_test_out, gld_test_test_out, 
                             nsd_test_out, sd_test_out, ymin=0.5, ymax=1):
    #order: 1)gld train 2) gld test 3) nsd 4) sd
    #return [prediction, allClasses_prediction_probabilities, accuracy, precision, recall, f1]
    #gather accuracies
    accuracies = [gld_train_test_out[2], gld_test_test_out[2],
                  nsd_test_out[2], sd_test_out[2]]
    #gather precisions
    precisions = [gld_train_test_out[3], gld_test_test_out[3],
                  nsd_test_out[3], sd_test_out[3]]
    #gather recalls
    recalls = [gld_train_test_out[4], gld_test_test_out[4],
                  nsd_test_out[4], sd_test_out[4]]
    #gather f1s
    f1s = [gld_train_test_out[5], gld_test_test_out[5],
                  nsd_test_out[5], sd_test_out[5]]
    #Dotplotting
    # Define custom x-axis labels
    x_labels = ["XLR_train", "XLR_test", "non-SD", "SD"]
    
    # Create x-values
    x_values = range(len(x_labels))
    
    # Original RGB color values
    acc_col = (0, 137, 255)
    f1_col = (0, 213, 255)
    prec_col = (255, 128, 0)
    rec_col = (255, 0, 0)
    
    # Convert to 0-1 range
    acc_col = (acc_col[0] / 255, acc_col[1] / 255, acc_col[2] / 255)
    f1_col = (f1_col[0] / 255, f1_col[1] / 255, f1_col[2] / 255)
    prec_col = (prec_col[0] / 255, prec_col[1] / 255, prec_col[2] / 255)
    rec_col = (rec_col[0] / 255, rec_col[1] / 255, rec_col[2] / 255)
    
    # Create a figure and axis
    fig, ax = plt.subplots()
    ax.set_ylim([ymin, ymax])
    
    # Plot the data series with connected lines    
    ax.plot(x_values, accuracies, marker='D', label='Accuracy', linestyle="solid", color=acc_col)
    ax.plot(x_values, f1s, marker='s', label='F1 Score', linestyle="dotted", color=f1_col)
    ax.plot(x_values, precisions, marker='o', label='Precision', linestyle="dashed", color=prec_col)
    ax.plot(x_values, recalls, marker='v', label='Recall', linestyle="dashdot", color=rec_col)
    
    # Set custom labels for the x-axis
    ax.set_xticks(x_values)
    ax.set_xticklabels(x_labels, **font) 

    # Set font for title, labels, and legend
    ax.set_xlabel('Dataset', **font)
    ax.set_ylabel('Value', **font)
    ax.legend(fontsize='large')
    
    # Add a legend to the plot
    ax.legend()
    
    figname = "{}_ConnectedScatter_Performance.pdf".format(modelname)
    plt.savefig(os.path.join(savepath, figname), dpi = 100)
    plt.close()
    return

def test_model(clf, model_output_path, disp_columns, training_columns,
               model_id, modelname, class_labels,
               x_train, x_test, x_nsd, x_sd,
               y_train, y_test, y_nsd, y_sd, 
               averaging, chrx_data_path):
    """
    Tests the model's performance and generates evaluation metrics, confusion matrices, and ROC curves.

    Parameters:
    clf: The trained model to be tested.
    model_output_path (str): The path to save the generated outputs.
    disp_columns: Features to display when plotting the mlp.
    training_columns: Features to build the ChrX datasets with predictions.
    model_id (str): The identifier of the model.
    modelname (str): The name of the model.
    class_labels (list): The list of class labels.
    x_train: The training input data.
    x_test: The testing input data.
    x_nsd: The noSegdup input data.
    x_sd: The SegDup input data.
    y_train: The training target labels.
    y_test: The testing target labels.
    y_nsd: The noSegdup target labels.
    y_sd: The SegDup target labels.
    averaging: Averaging method for multiclass metrics.

    Returns:
    The function does not return any value explicitly, but it generates evaluation outputs such as confusion matrices, ROC curves, and saves them in the specified model_output_path.
    """
    
    print("Making predicitions on XLR, noSegdup and SegDup datasets")
    #For tests on NSD and SD, predictions of classes -2 and 2+ will be merged into the classes -1 and +1 respectively
    #This is due to the lack of true -2 and 2+ exons in NSD and SD datasets
    
    class_labels_nsd_sd = [-1, 0, 1] #updating class labels
    
    gld_train, gld_test, nsd, sd = test_model_prediction_and_scores(clf, x_train, x_test,
                                                                    x_nsd, x_sd, 
                                                                    y_train, y_test, 
                                                                    y_nsd, y_sd, 
                                                                    avg = averaging)
   
    #building datasets with predictions and predictions probabilities
    train_with_predictions = pd.concat([pd.DataFrame(x_train, columns =training_columns ),
                                pd.Series(y_train, name='Class'),
                                pd.Series(gld_train[0], name = "{}_pred".format(model_id)),
                                pd.DataFrame(gld_train[1],
                                             columns = ["{}_pred_proba_{}".format(model_id, str(x)) for x in class_labels])], axis=1)
    test_with_predictions = pd.concat([pd.DataFrame(x_test, columns =training_columns ),
                                pd.Series(y_test, name='Class'),
                                pd.Series(gld_test[0], name = "{}_pred".format(model_id)),
                                pd.DataFrame(gld_test[1],
                                             columns = ["{}_pred_proba_{}".format(model_id, str(x)) for x in class_labels])], axis=1)
    
    nsd_with_predictions = pd.concat([pd.DataFrame(x_nsd, columns =training_columns ),
                                pd.Series(y_nsd, name='Class'),
                                pd.Series(nsd[0], name = "{}_pred".format(model_id)),
                                pd.DataFrame(nsd[1],
                                             columns = ["{}_pred_proba_{}".format(model_id, str(x)) for x in class_labels_nsd_sd])], axis=1)
    sd_with_predictions = pd.concat([pd.DataFrame(x_sd, columns =training_columns ),
                                pd.Series(y_sd, name='Class'),
                                pd.Series(sd[0], name = "{}_pred".format(model_id)),
                                pd.DataFrame(sd[1],
                                             columns = ["{}_pred_proba_{}".format(model_id, str(x)) for x in class_labels_nsd_sd])], axis=1)

    #get date of execution for prediction folder
    date = datetime.datetime.today().strftime('%Y,%m,%d')
    
    #creating folder for datasets with predictions
    data_dir = os.path.join(model_output_path, "ChrX data with prediction_{}".format(date))
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
        
    #saving dataframes with predictions and a date
    #TODO these datasets have their features scaled with the scaler
    #Will this be a problem? I certainly hope not but it will be
    train_with_predictions.to_csv(os.path.join(chrx_data_path,"XLR_train_pred.txt.gz"), compression = "gzip", sep = "\t", index = False)
    test_with_predictions.to_csv(os.path.join(chrx_data_path,"XLR_test_pred.txt.gz"), compression = "gzip", sep = "\t", index = False)
    nsd_with_predictions.to_csv(os.path.join(chrx_data_path,"NSD_pred.txt.gz"), compression = "gzip", sep = "\t", index = False)
    sd_with_predictions.to_csv(os.path.join(chrx_data_path,"SD_pred.txt.gz"), compression = "gzip", sep = "\t", index = False)
    
    # #if the model is an MLP, print it
    # if model_id.lower() in ["nn", "mlp"]:
    #     plot_nn(x_train, clf, model_output_path, disp_columns, modelname, gld_test[2]) #gld_test[2] = accuracy of model on test set
    
    
    #for each dataset do confusion matrix and roc curve
    
    #XLR train dataset
    cm_gld_train = confusion_matrix(y_train, gld_train[0], labels = class_labels)
    plot_confusion_matrix(modelname + "_XLR_train", cm_gld_train, class_labels, model_output_path, cmap=plt.cm.Blues)
    plot_roc_curve(modelname + "_XLR_train", model_output_path, y_train, gld_train[1])
    
    #XLR test dataset
    cm_gld_test = confusion_matrix(y_test, gld_test[0], labels = class_labels)
    plot_confusion_matrix(modelname + "_XLR_test", cm_gld_test, class_labels, model_output_path, cmap=plt.cm.Blues)
    plot_roc_curve(modelname + "_XLR_test", model_output_path, y_test, gld_test[1])
    
    #No segmental duplication dataset
    cm_nsd = confusion_matrix(y_nsd, nsd[0], labels = class_labels_nsd_sd)
    plot_confusion_matrix(modelname + "_NSD", cm_nsd, class_labels_nsd_sd, model_output_path, cmap=plt.cm.Blues)
    plot_roc_curve(modelname + "_NSD", model_output_path, y_nsd, nsd[1])
    
    #Segmental duplication dataset
    cm_sd = confusion_matrix(y_sd, sd[0], labels = class_labels_nsd_sd)
    plot_confusion_matrix(modelname + "_SD", cm_sd, class_labels_nsd_sd, model_output_path, cmap=plt.cm.Blues)
    plot_roc_curve(modelname + "_SD", model_output_path, y_sd, sd[1])
    
    #Calculate per-class metrics
    train_sc_metrics = single_class_metrics(train_with_predictions, class_labels, model_id)
    test_sc_metrics = single_class_metrics(test_with_predictions, class_labels, model_id)
    nsd_sc_metrics = single_class_metrics(nsd_with_predictions, class_labels_nsd_sd, model_id)
    sd_sc_metrics = single_class_metrics(sd_with_predictions, class_labels_nsd_sd, model_id)
    
    sc_metrics_all = {"XLR train": train_sc_metrics, 
                      "XLR test": test_sc_metrics,
                      "NSD": nsd_sc_metrics,
                      "SD": sd_sc_metrics}
    
    #metrics report (prints in terminal and saves on txt)
    extended_report(modelname, 
                    y_train, y_test, y_nsd, y_sd,
                    gld_train, gld_test, nsd, sd, model_output_path,
                    sc_metrics_all,
                    avg = averaging)
    #Dotplot for performance visualization on the four datasets
    plot_performance_dotplot(modelname, model_output_path, 
                             gld_train, gld_test,
                             nsd, sd)
    return