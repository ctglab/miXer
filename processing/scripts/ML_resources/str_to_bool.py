#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 12:24:57 2022

@author: elia
"""

# string to boolean converter function
import argparse

def str_to_bool(string, argname):
    """
    Converts a string representation of a boolean value to a boolean.

    Parameters:
    string (str): The string to be converted.
    argname (str): The name of the argument being checked.

    Returns:
    bool: The boolean value represented by the string.

    Raises:
    argparse.ArgumentTypeError: If the string does not represent a valid boolean value.
    """
    
    if isinstance(string, bool): #check if input is already a boolean
        return string
    else:
        lowercase_string = string.lower() #convert to lowercase to check for less words
        if lowercase_string in ["true", "yes", "y", "t", "whynot", "1", "ok"]:
            return True
        elif lowercase_string in ["false", "no", "n", "f", "nope", "0" "not"]:
            return False
        else:
            raise argparse.ArgumentTypeError("Boolean value expected for {}".format(argname))
  