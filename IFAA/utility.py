#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 19:39:15 2022

@author: ss
"""
import pandas as pd
import numpy as np
from itertools import compress

def colnames(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.columns.to_numpy(dtype = "U")

def ncol(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return len(x.columns)

def nrow(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return len(x.index)

def colSums(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis = 0)

def rowSums(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis = 1)

def r_in(x, y):
    x, y = np.array(x), np.array(y)
    return np.array([item in y for item in x])

def r_ni(x, y):
    x, y = np.array(x), np.array(y)
    return np.array([not item in y for item in x])

def rm(*args):
    del args

def slice_bool(x, y): 
    """Slicing list by a boolean list"""
    return list(compress(x, y))

def inv_bool(x):
    return [not i for i in x]

def cbind(x):
    return pd.concat(x, axis=1)

def which(x):
    return np.array([i for i, j in enumerate(x) if j])