#!/usr/bin/python
# -*- coding: utf-8 -*-
from loadData import *
import allUserDefinedFuncs
from metaData import *
import timeit
import numpy as np
import pandas as pd
import warnings



IFAA(haveDataM(), haveDataC(), linkIDname = "id")


def IFAA(
    MicrobData,
    CovData,
    linkIDname,
    testCov=None,
    ctrlCov=None,
    testMany=True,
    ctrlMany=False,
    nRef=40,
    nRefMaxForEsti=2,
    refTaxa=np.empty(0),
    adjust_method='BY',
    fdrRate=0.15,
    paraJobs=None,
    bootB=500,
    standardize=False,
    sequentialRun=False,
    refReadsThresh=0.2,
    taxkeepThresh=1,
    SDThresh=0.05,
    SDquantilThresh=0,
    balanceCut=0.2,
    seed=1,
    ):
    
    testCov = np.array(testCov)
    ctrlCov = np.array(ctrlCov)
    linkIDname = np.array(linkIDname)

    allFunc = allUserDefinedFuncs.allUserFunc()
    results = {}
    
    start = timeit.default_timer()
    stop = timeit.default_timer()
    print('Time: ', stop - start) 
    runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                     linkIDname=linkIDname,taxkeepThresh=taxkeepThresh,
                   testCov=testCov,
                   ctrlCov=ctrlCov,testMany=testMany,
                   ctrlMany=ctrlMany)
    
    data=runMeta['data']
    results['covariatesData']=runMeta['covariatesData']
    binaryInd=runMeta['binaryInd']
    covsPrefix=runMeta['covsPrefix']
    Mprefix=runMeta['Mprefix']
    testCovInd=runMeta['testCovInd']
    testCovInOrder=runMeta['testCovInOrder']
    testCovInNewNam=runMeta['testCovInNewNam']
    ctrlCov=runMeta['ctrlCov']
    microbName=runMeta['microbName']
    newMicrobNames=runMeta['newMicrobNames']
    results['covriateNames']=runMeta['xNames']
    del(runMeta)
    
    if (refTaxa is not None) and (len(refTaxa)>0):
        if sum(r_in(refTaxa, microbName)) != len(refTaxa):
            raise Exception("""
                             Error: One or more of the specified reference taxa in phase 1 have no sequencing reads 
                             or are not in the data set. Double check the names of the reference taxa and their 
                             sparsity levels.""")
    if nRefMaxForEsti<2:
        nRefMaxForEsti = 2
        warnings.warn("Warning: Needs at least two final reference taxon for estimation.")

    if nRef>len(microbName):
        raise Exception("Error: number of random reference taxa can not be larger than the total number of taxa in the data. Try lower nRef")
        
    refTaxa_newNam=newMicrobNames[r_in(microbName, refTaxa)]

    results['analysisResults']

