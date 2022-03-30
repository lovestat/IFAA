#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 22:07:38 2022

@author: ss
"""
## Original screen function to select independent reference taxa

import numpy as np
from utility import *
import timeit
import multiprocessing as mp
import math
import joblib
from dataInfo import *
from dataRecovTrans import *
from AIcalcu import *

from sklearn.linear_model import LassoCV

from statsmodels.api import OLS
from statsmodels.stats.multitest import multipletests


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
  maxDimensionScr=434*5*10**5):
    
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
    
    countOfSelec=np.zeros(nAlphaSelec)
    
    # overwrite nRef if the reference taxon is specified
    nRef=len(refTaxa)
    
    startT1=timeit.default_timer()
    if paraJobs is None:
        availCores = mp.cpu_count()
        if isinstance(availCores, int): 
            paraJobs = max(1, availCores-2)
    
    batch=paraJobs
    forLoopN=math.ceil(nRef/batch)
    
    if not sequentialRun:
        print(paraJobs, " parallel jobs are registered for analyzing ", nRef, " reference taxa in Phase 1")

        
i = 1

def foreachUnitRun(i):
    
    ii=which(taxaNames==refTaxa[i])
    dataForEst=dataRecovTrans(data=data,
                              ref=refTaxa[i],
                              Mprefix=Mprefix,
                              covsPrefix=covsPrefix)

    xTildLongTild_i=dataForEst['xTildalong']
    yTildLongTild_i=dataForEst['UtildaLong']
    rm(dataForEst)
    
    maxSubSamplSiz=np.min( (50000.0, np.floor(maxDimensionScr / xTildLongTild_i.shape[1]))).astype(int)
    
    nToSamplFrom = xTildLongTild_i.shape[0]

    subSamplK=np.ceil(nToSamplFrom/maxSubSamplSiz).astype(int)
    
    if subSamplK==1 : maxSubSamplSiz=nToSamplFrom

    nRuns=np.ceil(subSamplK/3).astype(int)
    
    for k in range(nRuns):
        rowToKeep = np.random.choice(nToSamplFrom, maxSubSamplSiz, replace = False)
        
        x = xTildLongTild_i[rowToKeep, :]
        y = yTildLongTild_i[rowToKeep]
        
        if x.shape[0] > (3 * x.shape[1]):
            Penal.i=runlinear(x=x,y=y, nPredics=nPredics)
            BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
            EstNoInt.k<-abs(Penal.i$coef_est_noint)
        else:
            Penal.i=runGlmnet(x=x,y=y, nPredics=nPredics, standardize=standardize)
            BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
            EstNoInt.k<-abs(Penal.i$betaNoInt)
        
        
def runlinear(
        x,
        y,
        nPredics,
        fwerRate=0.25,
        adjust_method="fdr"):
    
    results={}
    nBeta=x.shape[1]
    nObsAll=len(y)
    print("length of y: ",len(y))
    
    #lm_model = LinearRegression(fit_intercept = False)
    #lm_res = lm_model.fit(x, y)
    #lm_res.coef_
    
    lm_res = OLS(y, x).fit()
    p_value_est = lm_res.pvalues
    disc_index = np.arange(0, len(p_value_est), (nPredics+1) )
    p_value_est_noint = np.delete(p_value_est, disc_index, axis=0)
    
    ## this method automatically convert over 1 values to 1
    p_value_est_noint_adj = multipletests(p_value_est_noint, 
                                          alpha = 0.05,
                                          method = adjust_method)[1]

    coef_est = np.abs(lm_res.params)
    disc_index = np.arange(0, len(p_value_est), (nPredics+1) )
    ## NA coef is snot considered here
    coef_est_noint = np.delete(coef_est, disc_index, axis=0)
    
    # return
    results['betaNoInt']=p_value_est_noint_adj<fdrRate
    results['betaInt']=p_value_est
    results['coef_est_noint']=coef_est_noint


    return results
    

  standardize=False
  family="gaussian"
  nfolds=10
  lambda_min_ratio=0.05
  nLam=100
  intercept=True
  zeroSDCut=10**(-20)
  
  

def runGlmnet(
  x,
  y,
  nPredics,
  standardize=False,
  family="gaussian",
  nfolds=10,
  lambda_min_ratio=0.05,
  nLam=100,
  intercept=True,
  zeroSDCut=10**(-20)
):
    results={}
    nBeta=x.shape[1]
    nObsAll=len(y)
    
    # remove near constant x columns
    sdX=np.std(x, axis = 0)
    xWithNearZeroSd=which(sdX<=zeroSDCut)
    if len(xWithNearZeroSd)>0 :
        x = np.delete(x, xWithNearZeroSd, axis=1)
    rm(sdX)
    
    # calculate lambda max
    if family=="gaussian" :
        lamMax = np.max( np.abs( (x*y[:, np.newaxis]).sum(0) ) ) / nObsAll        
        lamVec = np.linspace(lamMax, 0, nLam+1)[0:nLam]  
        cvStartTime = timeit.default_timer()
        
        cvStartTimei = timeit.default_timer()
        
        
        
        cvResul=LassoCV(alphas=lamVec, 
                        fit_intercept=intercept,
                        normalize = standardize,
                        cv = nfolds, 
                        n_jobs = -1).fit(x.copy(), y.copy())
        cvResul.fit(x, y)
        
        cvExeTimei= (timeit.default_timer() - cvStartTimei)/60

        cvResul.alpha
        lamOpi=as.numeric(cvResul$lambda.min)
        cvm=as.vector(cvResul$cvm)*nObsAll
        nLamUsed=length(as.vector(cvResul$lambda))

      rm(cvResul)
    
    