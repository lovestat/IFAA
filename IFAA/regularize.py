#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:49:50 2022

@author: jin
"""

map(eval, """data=data,Mprefix=Mprefix,
                 covsPrefix=covsPrefix,
                 binPredInd=binaryInd""".split(","))



def Regulariz(
        data,
        testCovInd,
        testCovInOrder,
        testCovInNewNam,
        microbName,
        nRef,
        nRefMaxForEsti,
        refTaxa,
        paraJobs,
        binaryInd,
        covsPrefix,
        Mprefix,
        fwerRate,
        bootB,
        standardize,
        sequentialRun,
        allFunc=allFunc,
        refReadsThresh,
        SDThresh,
        SDquantilThresh,
        balanceCut,
        adjust_method,
        seed):
    results = {}
    regul_start_time = timeit.default_timer()
    
    nTestCov=len(testCovInd)
    dataSparsCheck(data=data,Mprefix=Mprefix)

    # load abundance data info
    
    binCheck=data.loc[:, testCovInNewNam].nunique()
    binaryInd<-which(binCheck==2)
    
    data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)
    
    
    
    
    