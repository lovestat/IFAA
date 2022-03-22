# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:16:52 2022

@author: ss
"""
import pandas as pd
import warnings


def colnames(x):
    if not isinstance(x, pd.core.frame.DataFrame): 
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.columns.tolist()

def r_in(x, y):
    return pd.Series(x).isin(y)

pd.Series(testCov + ctrlCov).isin(CovData.columns.tolist())



def metaData(
        MicrobData,
        CovData,
        linkIDname,
        taxkeepThresh,
        testCov=None,
        ctrlCov=None,                  
        testMany=True,
        ctrlMany=True,
        MZILN=True):
    
    results={}
    
    if not linkIDname:
        raise Exception("linkIDname is missing.")
        
    if ( len(testCov)>0 and len(ctrlCov)>0 ):
        if sum(r_in(testCov + ctrlCov, colnames(CovData))) != len(testCov + ctrlCov):
            raise Exception("Error: some covariates are not available in the data.")
  
    if set(testCov).intersection(ctrlCov) :
        warnings.warn('Variables appeared in both testCov list and ctrlCov list will be treated as testCov.')
        
    # read microbiome data
    MdataWithId=MicrobData
    if (MdataWithId[MdataWithId.columns.drop(linkIDname)].values < 0).any():
        raise Exception("Microbiome data contains negative values.")
    
    if MdataWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn("There are over 80% missing values for the linkId variable in the Microbiome data file. Double check the data format.")
    
    # read covariate data
    CovarWithId=CovData
    if CovarWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn("There are over 80% missing values for the linkId variable in the covariates data file. Double check the data format.")
    
    Covariates1 = CovarWithId[CovarWithId.columns.drop(linkIDname)]
    
    # determine testCov and ctrlCov
    if len(testCov) == 0:
        if not testMany:
            raise Exception("No covariates are specified for estimating associations of interest.")
        else:
            print('Associations are being estimated for all covariates since no covariates are specified for testCov.')
            testCov = Covariates1.columns.tolist()
    results['testCov'] = testCov
    
    xNames=Covariates1.columns.tolist()
    del Covariates1
    
    if len(ctrlCov) == 0 and ctrlMany:
        print('No control covariates are specified, all variables except testCov are considered as control covariates.')
        ctrlCov = set(xNames).difference(testCov)
    ctrlCov = list(ctrlCov.difference(testCov))
    results['ctrlCov'] = ctrlCov
    
    # merge data to remove missing
    CovarWithId1=CovarWithId[ [linkIDname]+testCov+ctrlCov ]
    allRawData = CovarWithId1.merge(MdataWithId, on = linkIDname).dropna()
    CovarWithId=allRawData[,(colnames(allRawData)%in%colnames(CovarWithId1)),drop=FALSE]
    Covariates=CovarWithId[,!colnames(CovarWithId)%in%linkIDname,drop=FALSE]
    

metaData(haveDataM(), haveDataC(), "id", None, ["v1", "v2"], ["v3"])
linkIDname = "id"
testCov = ["v1", "v2"]
ctrlCov = ["v3"]

set(testCov)

