#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 12:23:00 2022

@author: elia
"""

#function that splits lists in chunks
#used for multithreading

def split(list_a, chunk_size):
    """
    Splits a list into smaller sublists of a specified chunk size.

    Parameters:
    list_a (list): The list to be split.
    chunk_size (int): The size of each sublist.

    Yields:
    list: A sublist of the original list with the specified chunk size.
    """

    for i in range(0, len(list_a), chunk_size):
        yield list_a[i:i + chunk_size]