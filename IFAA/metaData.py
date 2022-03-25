# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 12:16:52 2022

@author: ss
"""
import pandas as pd
import numpy as np
import warnings
from utility import *



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
        if sum(r_in(np.concatenate((testCov, ctrlCov)), colnames(CovData))) != len(np.concatenate((testCov, ctrlCov))):
            raise Exception("Error: some covariates are not available in the data.")
  
    if sum( r_in(testCov, ctrlCov) ) > 0:
        warnings.warn('Variables appeared in both testCov list and ctrlCov list will be treated as testCov.')
        
    # read microbiome data
    MdataWithId=MicrobData
    if( len(colnames(MdataWithId)) != ncol(MdataWithId) ): 
        raise Exception("Microbiome data lack variable names.")
    
    if (MdataWithId.loc[:, MdataWithId.columns != linkIDname].values < 0).any():
        raise Exception("Microbiome data contains negative values.")
    
    if MdataWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn("There are over 80% missing values for the linkId variable in the Microbiome data file. Double check the data format.")
    
    
    
    # read covariate data
    CovarWithId=CovData
    if CovarWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn("There are over 80% missing values for the linkId variable in the covariates data file. Double check the data format.")
    
    Covariates1 = CovarWithId.loc[:, CovarWithId.columns != linkIDname]
    
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
        ctrlCov= xNames[r_ni(xNames, testCov)]  
    ctrlCov = ctrlCov[r_ni(ctrlCov, testCov)]
    results['ctrlCov'] = ctrlCov
    
    # merge data to remove missing
    CovarWithId1=CovarWithId[ np.hstack([linkIDname, testCov, ctrlCov]) ]
    allRawData = CovarWithId1.merge(MdataWithId, on = linkIDname.tolist()).dropna()
    CovarWithId = allRawData.loc[:, r_in(colnames(allRawData), colnames(CovarWithId1))]
    Covariates=CovarWithId.loc[:, r_ni(colnames(CovarWithId), [linkIDname])]
    rm(CovarWithId1)
    
    if not all(Covariates.apply(lambda a: a.dtype.kind in 'if')):
        raise Exception("There are non-numeric variables in the covariates for association test.")
        
    MdataWithId=allRawData.loc[:, r_in(colnames(allRawData), colnames(MdataWithId)) ]
    Mdata_raw=MdataWithId.loc[:, r_ni(colnames(MdataWithId), linkIDname)]
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
    newMicrobNames1=np.array(["microb" + str(i+1) for i in range(len(microbName))])
    newMicrobNames=newMicrobNames1
    results['Mprefix']="microb"    
    Mdata=Mdata.rename(columns=dict( zip(microbName1, newMicrobNames) ))
    
    MdataWithId_new=cbind([MdataWithId.loc[:,linkIDname], Mdata])
    
    results['microbName']=microbName
    results['newMicrobNames']=newMicrobNames
    
    if Covariates.isna().sum().sum() >0 :
        print("Samples with missing covariate values are removed from the analysis.")
    
    ## to add non-numeric add
    
    xNames=colnames(Covariates)
    nCov=len(xNames)
    
    ## to add binary check
    
    results['nBinVars']=0
    binaryInd=None
    results['varNamForBin']=None
    
    results['binaryInd']=None
    results['xNames']=colnames(Covariates)
    xNewNames=np.array(["x" + str(i+1) for i in range(len(xNames))])
    Covariates = Covariates.rename(columns=dict( zip(colnames(Covariates), xNewNames) ))
    results['covsPrefix']="x"
    results['xNewNames']=xNewNames
    
    results['testCovInd']=which(r_in(results['xNames'],testCov))
    
    results['testCovInOrder']=results['xNames'][results['testCovInd']]
    results['testCovInNewNam']=results['xNewNames'][results['testCovInd']]
    del(xNames,xNewNames)
      
    CovarWithId_new=cbind([CovarWithId.loc[:,linkIDname], Covariates])
    data = MdataWithId_new.merge(CovarWithId_new, on = linkIDname.tolist())
    dataOmit = data.dropna()
    
    
    results['covariatesData']=CovarWithId_new
    results['covariatesData'].rename(columns=dict( zip(results['covariatesData'], np.insert(results['xNames'], 0, linkIDname) ) ))
    del(MdataWithId_new,CovarWithId_new)

    Mdata_omit=dataOmit.loc[:,np.array(newMicrobNames1)]
    
    # check taxa with zero or 1 read again after all missing data removed
    # to add
    
    results['data']=dataOmit
    
    return results
    
    
    
    
