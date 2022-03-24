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
    return x.columns.to_numpy()

def ncol(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return len(x.columns)

def colSums(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis = 0)

def rowSums(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis = 1)

def r_in(x, y):
    return pd.Series(x).isin(y).tolist()

def r_ni(x, y):
    return [not item for item in r_in(x, y)]

def rm(x):
    del x

def slice_bool(x, y): 
    """Slicing list by a boolean list"""
    return list(compress(x, y))

def inv_bool(x):
    return [not i for i in x]

def cbind(x):
    return pd.concat(x, axis=1)

def which(x):
    return np.array([i for i, j in enumerate(x) if j])