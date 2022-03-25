#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 19:37:49 2022

@author: jin
"""

def dataInfo(
    data,
    Mprefix,
    covsPrefix,
    binPredInd,
    refReadsThresh,
    SDThresh,
    SDquantilThresh,
    balanceCut,
    qualifyRefTax=False,
):
    results={}
    
    # get the original sample size
    nSub=nrow(data)
    MVarNamLength=len(Mprefix)
    
    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    nTaxa = len(which(microPositions))
    taxaNames = data.columns[microPositions]
    rm(microPositions)
    
    ## to add if qualifyreftax 
    
    # get predictor data
    predNames=data.columns[data.columns.str.startswith(covsPrefix)]
    nPredics=len(predNames)
    
    ## to add if qualifyreftax 
    
    # return
    results['taxaNames']=taxaNames
    rm(taxaNames)
    results['predNames']=predNames
    rm(predNames)
    results['nTaxa']=nTaxa
    results['nSub']=nSub
    results['nPredics']=nPredics
    return results
    
    
    
    
