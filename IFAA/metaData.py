# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:16:52 2022

@author: ss
"""
import pandas as pd
import warnings
from utility import *

linkIDname = "id"
testCov = ["v1", "v2"]
ctrlCov = ["v3"]
taxkeepThresh = 0

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
  
    if sum( r_in(testCov, ctrlCov) ) > 0:
        warnings.warn('Variables appeared in both testCov list and ctrlCov list will be treated as testCov.')
        
    # read microbiome data
    MdataWithId=MicrobData
    if( len(colnames(MdataWithId)) != ncol(MdataWithId) ): 
        raise Exception("Microbiome data lack variable names.")
    
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
            testCov = colnames(Covariates1)
    results['testCov'] = testCov
    
    xNames=colnames(Covariates1)
    rm(Covariates1)
    
    if len(ctrlCov) == 0 and ctrlMany:
        print('No control covariates are specified, all variables except testCov are considered as control covariates.')
        ctrlCov= slice_bool(xNames, r_ni(xNames, testCov))   
    ctrlCov = slice_bool(ctrlCov, r_ni(ctrlCov, testCov))
    results['ctrlCov'] = ctrlCov
    
    # merge data to remove missing
    CovarWithId1=CovarWithId[ [linkIDname]+testCov+ctrlCov ]
    allRawData = CovarWithId1.merge(MdataWithId, on = linkIDname).dropna()
    CovarWithId = allRawData.loc[:, r_in(colnames(allRawData), colnames(CovarWithId1))]
    Covariates=CovarWithId.loc[:, r_ni(colnames(CovarWithId), [linkIDname])]
    rm(CovarWithId1)
    
    if not all(Covariates.apply(lambda a: a.dtype.kind in 'if')):
        raise Exception("There are non-numeric variables in the covariates for association test.")
        
    MdataWithId=allRawData.loc[:, r_in(colnames(allRawData), colnames(MdataWithId)) ]
    Mdata_raw=MdataWithId.loc[:, r_ni(colnames(MdataWithId), [linkIDname])]
    rm(allRawData)
    
    # check zero taxa and subjects with zero taxa reads
    numTaxaNoReads = sum(colSums(Mdata_raw != 0) <= taxkeepThresh)
    if numTaxaNoReads>0:
        Mdata_raw=Mdata_raw.loc[:,colSums(Mdata_raw != 0) > taxkeepThresh]
        print("There are ",numTaxaNoReads," taxa without any sequencing reads before data merging, and excluded from the analysis")
    rm(numTaxaNoReads)
    
    numSubNoReads=sum(rowSums(Mdata_raw!=0)<=1)
    if numSubNoReads>0:
        print("There are ",numSubNoReads," subjects with zero or one sequencing read and excluded from the analysis")
        subKeep= inv_bool(rowSums(Mdata_raw!=0)<=1)
        Mdata_raw=Mdata_raw.loc[subKeep, :]
        MdataWithId=MdataWithId.loc[subKeep, :]
        rm(subKeep)
    rm(numSubNoReads)
    
    Mdata=Mdata_raw
    rm(Mdata_raw)
    
    microbName1=colnames(Mdata)
    microbName=microbName1
    newMicrobNames1=["microb" + str(i+1) for i in range(len(microbName))]
    newMicrobNames=newMicrobNames1
    results['Mprefix']="microb"    
    Mdata=Mdata.rename(columns=dict( zip(microbName1, newMicrobNames) ))
    
    MdataWithId_new=cbind([MdataWithId.loc[:,linkIDname], Mdata])
    
    results['microbName']=microbName
    results['newMicrobNames']=newMicrobNames
    
    if Covariates.isna().sum().sum() >0 :
        print("Samples with missing covariate values are removed from the analysis.")
    
 

    
    # dictionary
Data = {'Name': ['GeeksForGeeks','Python'],
          'Unique ID': ['900','450']}
 
# create a dataframe object
df = pd.DataFrame(Data)

pd.to_numeric()
    
df.
    
    
    
    
    