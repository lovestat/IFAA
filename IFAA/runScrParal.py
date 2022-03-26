#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 21:32:12 2022

@author: jin
"""
import numpy as np

def runScrParal(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  refTaxa,
  maxDimensionScr=0.8*434*10*10**4,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
  adjust_method,
  seed):
    
    # load data info
    basicInfo=dataInfo(data=data,
                       Mprefix=Mprefix,
                       covsPrefix=covsPrefix,
                       binPredInd=binPredInd,
                       refReadsThresh=refReadsThresh,
                       SDThresh=SDThresh,
                       SDquantilThresh=SDquantilThresh,
                       balanceCut=balanceCut,
                       qualifyRefTax=True)
    
    taxaNames=basicInfo['taxaNames']
    nTaxa=basicInfo['nTaxa']
    nPredics=basicInfo['nPredics']
    nSub=basicInfo['nSub']
    predNames=basicInfo['predNames']
    
    results['goodRefTaxaCandi']=basicInfo['goodRefTaxaCandi']
    rm(basicInfo)
    
    nNorm=nTaxa-1
    nAlphaNoInt=nPredics*nNorm
    nAlphaSelec=nPredics*nTaxa
    
    # make reference taxa list
    if len(refTaxa)<nRef:
        if len(seed)>0:
            np.random.seed(seed)
            
        taxon_to_be_sample = results['goodRefTaxaCandi'][r_ni(results['goodRefTaxaCandi'], refTaxa)]
        num_to_be_sample = (nRef-len(refTaxa))
        
        if num_to_be_sample >= len(taxon_to_be_sample):
            num_to_be_sample = len(taxon_to_be_sample)
        print("The number of candidate reference taxon is smaller than the number of taxon required in phase 1. The number of taxon was set to be ",num_to_be_sample)
        
        refTaxa_extra = np.random.choice(taxon_to_be_sample,num_to_be_sample, replace=False)
        refTaxa = np.hstack( (refTaxa,refTaxa_extra))
        results['refTaxa'] = refTaxa
    
        if len(refTaxa)==0:
            raise Exception("No candidate reference taxon is available. Please try to lower the reference taxon boundary.")
            
        
    if len(refTaxa) >= nRef:
        if len(seed) > 0:
            np.random.seed(seed)
        refTaxa=np.random.choice(refTaxa, nRef, replace=True)
        results['refTaxa'] = refTaxa
    
    ## run original data screen
    screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,
                           paraJobs=paraJobs,
                           allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           adjust_method=adjust_method,
                           seed=seed)
    
    
