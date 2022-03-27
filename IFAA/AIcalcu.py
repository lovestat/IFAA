#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 00:06:30 2022

@author: jin
"""

import numpy as np

def AIcalcu(
        data,
        ref,
        Mprefix,
        covsPrefix):
    results = {}
    
    # get the original sample size
    nSub=nrow(data)
    MVarNamLength=len(Mprefix)
    
    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    nTaxa = len(which(microPositions))
    taxaNames = data.columns[microPositions]
    rm(microPositions)
    
    # rearrange taxa names
    otherTaxaNames=taxaNames[r_ni(taxaNames, ref)]
    taxaNames=np.hstack([otherTaxaNames,ref])

    # get predictor data
    predNames=data.columns[data.columns.str.startswith(covsPrefix)].to_numpy()
    nPredics=len(predNames)

    # taxa data
    w=data.loc[:,taxaNames]
    
    # extract x data
    xData=data.loc[:,predNames]
    rm(data,predNames)

    # transform data using log-ratio, creat Ai and Li
    l=np.empty(nSub)
    lLast=np.empty(nSub)
    taxa_non0={}
    taxa_0={}
    logRatiow={}
    A={}
    
    for i in range(nSub):
        taxa_nonzero=which(w.iloc[i,:]!=0)
        lLast[i]=np.max(taxa_nonzero)
        taxa_zero=which(w.iloc[i,:]==0)
        taxa_non0[i]=w.iloc[i, taxa_nonzero]
        taxa_0[i]=w.iloc[i,taxa_zero]
        if len(taxa_nonzero) > 0:
            last_nonzero=np.max(taxa_nonzero)
            logwi=np.log(w.iloc[i,taxa_nonzero])
            l[i]=len(logwi)
            if l[i]>1 :
                logRatiow[i]=logwi[:-1:]-logwi[-1]
                zero_m=np.zeros( (int(l[i])-1, nNorm) )
                if last_nonzero == nTaxa:
                    aRow = np.arange(int(l[i])-1)
                    aCol = taxa_nonzero[:-1:]
                    zero_m[aRow, aCol]=1
                else:
                    aRow = np.arange(int(l[i])-1)
                    aCol = taxa_nonzero[:-1:]
                    zero_m[aRow, aCol]=1
                    zero_m[:, int(taxa_nonzero[int(l[i])-1])-1 ] = -1
                    A[i]=zero_m
                    rm(zero_m)
            else:
                logRatiow[i]=None
                A[i]=None
        else:
            l[i]=0
            logRatiow[i]=None
            A[i]=None
    
    # obtain the list of samples whose have at least 2 non-zero taxa
    twoList=which(l>1)
    lengthTwoList=len(twoList)
    
    rm(w)
    
    results['xData']=xData
    rm(xData)
    
    results['logRatiow']=logRatiow
    rm(logRatiow)
    results['A']=A
    rm(A)
    results['twoList']=twoList
    rm(twoList)
    results['taxaNames']=taxaNames
    rm(taxaNames)
    results['lengthTwoList']=lengthTwoList
    results['lLast']=lLast
    results['l']=l
    results['nTaxa']=nTaxa
    results['nNorm']=nNorm
    results['nSub']=nSub
    results['nPredics']=nPredics
    return results
        
        