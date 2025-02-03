#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 09:37:11 2022

@author: elia
"""

#function that plots the histogram of a dataframe column
import os
import matplotlib
import matplotlib.pyplot as plt
#setting matplotlib backend to work in non-interactive mode
matplotlib.use('Agg')
plt.rcParams.update({'font.size': 15})

def plot_df_hist(df, column_to_plot_name, df_name, output_path):
    """
    Plots a histogram for a specific column in a DataFrame.

    Args:
        df (pandas.DataFrame): Input DataFrame.
        column_to_plot_name (str): Name of the column to plot.
        df_name (str): Name of the DataFrame.
        output_path (str): Output directory path for saving the histogram plot.

    Returns:
        None
    """
    
    figname = df_name + column_to_plot_name + "_histogram.png"
    fig, ax = plt.subplots(figsize = (10,5))
    print("{} max: {} | Min: {}".format(column_to_plot_name, round(df[column_to_plot_name].max(), 2), round(df[column_to_plot_name].min(), 2)))
    ax = df[column_to_plot_name].hist(bins = 200)
    fig.savefig(os.path.join(output_path, figname), dpi = 300)
    return