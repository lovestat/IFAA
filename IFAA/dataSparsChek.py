#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 18:10:28 2022

@author: jin
"""

from utility import *

def dataSparsCheck(data, Mprefix):
    results = {}
    
    # get the original sample size
    nSub=nrow(data)
    MVarNamLength=len(Mprefix)
    
    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    taxaNames = data.columns[microPositions]
    rm(microPositions)
    
    w=data.loc[:, taxaNames]
    rm(data,taxaNames)
    overallSparsity= np.round(100*np.mean(w.values == 0), 2)
    print(overallSparsity, "percent of microbiome sequencing reads are zero")

    # check zero taxa and subjects with zero taxa reads
    numTaxaNoReads=sum(colSums(w)==0)
    if numTaxaNoReads>0:
        print("There are ",numTaxaNoReads," taxa without any sequencing reads and excluded from the analysis")
    rm(numTaxaNoReads)
    
    numSubNoReads=sum(rowSums(w)==0)
    if numSubNoReads>0:
        print("There are ",numSubNoReads," subjects without any sequencing reads and excluded from the analysis.")
    rm(numSubNoReads,w)
    