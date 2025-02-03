#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 17:29:49 2022

@author: elia
"""

def get_mean(l):
    """
    Calculates the mean of a given list of numbers.

    Args:
        l (list): List of numbers.

    Returns:
        float: Mean of the numbers in the list.
    """
    return sum(l) / len(l)


def get_var(l):
    """
    Calculates the variance of a given list of numbers.

    Args:
        l (list): List of numbers.

    Returns:
        float: Variance of the numbers in the list.
    """
    mean = sum(l) / len(l)
    return sum((x - mean) ** 2 for x in l) / len(l)



def get_std(l):
    """
    Calculates the standard deviation of a given list of numbers.

    Args:
        l (list): List of numbers.

    Returns:
        float: Standard deviation of the numbers in the list.
    """
    variance = get_var(l)
    return variance ** 0.5
