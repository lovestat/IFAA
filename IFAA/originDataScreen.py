#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 22:07:38 2022

@author: ss
"""
## Original screen function to select independent reference taxa

import numpy as np
from utility import *

def originDataScreen(
  data,
  testCovInd,
  nRef,
  paraJobs,
  refTaxa,
  standardize,
  sequentialRun,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  adjust_method,
  seed,
  maxDimensionScr=434*5*10^5):
    
    results = {}
    
    # load data info
    basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                   covsPrefix=covsPrefix,
                   binPredInd=binPredInd)
    
    taxaNames=basicInfo['taxaNames']
    nTaxa=basicInfo['nTaxa']
    nPredics=basicInfo['nPredics']
    rm(basicInfo)
    
    nNorm=nTaxa-1
    nAlphaNoInt=nPredics*nNorm
    nAlphaSelec=nPredics*nTaxa
    
    countOfSelec=rep(0, nAlphaSelec)
    
    # overwrite nRef if the reference taxon is specified
    nRef=len(refTaxa)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
