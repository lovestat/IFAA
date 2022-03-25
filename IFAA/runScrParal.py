#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 21:32:12 2022

@author: jin
"""

def runScrParal(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  refTaxa,
  maxDimensionScr=0.8*434*10*10^4,
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
