#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 00:06:08 2022

@author: jin
"""
import numpy as np
from AIcalcu import *

def dataRecovTrans(
  data,
  ref,
  Mprefix,
  covsPrefix,
  xOnly=False,
  yOnly=False
):
    results={}

    # load A and log-ratio transformed RA
    data_and_init=AIcalcu(data=data,ref=ref,Mprefix=Mprefix,
                          covsPrefix=covsPrefix)
    rm(data)

    taxaNames=data_and_init['taxaNames']
    A=data_and_init['A']
    logRatiow=data_and_init['logRatiow']
    nSub=data_and_init['nSub']
    nTaxa=data_and_init['nTaxa']
    xData=data_and_init['xData']
    nPredics=data_and_init['nPredics']
    twoList=data_and_init['twoList']
    lLast=data_and_init['lLast']
    L=data_and_init['l']
    lengthTwoList=data_and_init['lengthTwoList']
    rm(data_and_init)
    
    nNorm=nTaxa-1
    xDimension=nPredics+1 # predictors+intercept
    
    # create omegaRoot
    omegaRoot={}
    for j in range(lengthTwoList):
        i=twoList[j]
        if lLast[i]==nTaxa:
            omegaRoot[i]=np.eye( int(L[i]-1) )
        else: 
            if L[i]==2:
                omegaRoot[i]=np.sqrt(0.5) 
            else:
                dim=L[i]-1
                a=(1+(dim-2)/2)/(1/2*(1+(dim-1)/2))
                b=-1/(1+(dim-1)/2)
            
                # calculate the square root of omega assuming it is exchangeable
                aStar=dim**2/((dim-1)**2)
                bStar=b*(dim-2)/(dim-1)-a*(dim**2-2*dim+2)/((dim-1)**2)
                cStar=(0.5*b-0.5*a*(dim-2)/(dim-1))**2
                cSquare=(-bStar+np.sqrt(bStar**2-4*aStar*cStar))/(2*aStar)
                
                if cSquare<0 :
                    raise Exception("no solution for square root of omega")
                d=np.sqrt((0.5*a-cSquare)/(dim-1))
                
                if  d is None :
                    raise Exception("no solution for off-diagnal elements for square root of omega")
                c=(0.5*b-(dim-2)*(d**2))/(2*d)
                omegaRoot[i]=-((c-d)*np.eye(int(dim))+d*np.ones( (int(dim), int(dim)) ))
     
        rm(L,lLast)
    
    if xOnly :
        # create X_i in the regression equation using Kronecker product
        xDataWithInter= xData.copy()   
        xDataWithInter.insert(0, "Intercept", 1)
        xDataWithInter=xDataWithInter.to_numpy()
        rm(xData)
        
        for j in range(lengthTwoList):
            i=twoList[j]
            xInRegres_i=np.kron(np.eye(nNorm),xDataWithInter[i,:])
            xDataTilda_i=omegaRoot[i] @ A[i] @ xInRegres_i
            rm(xInRegres_i)
            
            if j == 0:
                xTildalong=xDataTilda_i
            else:
                xTildalong=np.vstack((xTildalong,xDataTilda_i))
        
        rm(xDataWithInter,xDataTilda_i,omegaRoot,logRatiow)

        results['xTildalong']=xTildalong
        rm(xTildalong)
        return results
    
    if yOnly:
        for j in range(lengthTwoList):
            i=twoList[j]
            Utilda_i=omegaRoot[i] @ logRatiow[i]
            if j == 0:
                UtildaLong = Utilda_i
            else:
                UtildaLong=np.hstack( (UtildaLong,Utilda_i) )
    
        rm(omegaRoot,logRatiow,Utilda.i)
        results['UtildaLong']=UtildaLong
        rm(UtildaLong)
        return results
    
    # create X_i in the regression equation using Kronecker product
    xDataWithInter = xData.copy()   
    rm(xData)
    xDataWithInter.insert(0, "Inter", 1)
    xDataWithInter = xDataWithInter.to_numpy()
    
    for j in range(lengthTwoList):
        i=twoList[j]
        xInRegres_i=np.kron(np.eye(nNorm),xDataWithInter[i,:])
        xDataTilda_i=omegaRoot[i] @ A[i] @ xInRegres_i
        rm(xInRegres_i)
        
        if j == 0:
            xTildalong=xDataTilda_i
        else:
            xTildalong=np.vstack((xTildalong,xDataTilda_i))
    
    rm(xDataWithInter,xDataTilda_i)
    
    for j in range(lengthTwoList):
        i=twoList[j]
        Utilda_i=omegaRoot[i] @ logRatiow[i]
        if j == 0:
            UtildaLong = Utilda_i
        else:
            UtildaLong=np.hstack( (UtildaLong,Utilda_i) )

    rm(omegaRoot,logRatiow,Utilda_i)
    
    # return objects
    results['UtildaLong']=UtildaLong
    rm(UtildaLong)
    results['xTildalong']=xTildalong
    rm(xTildalong)
    results['taxaNames']=taxaNames
    rm(taxaNames)
    return results
