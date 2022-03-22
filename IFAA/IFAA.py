#!/usr/bin/python
# -*- coding: utf-8 -*-
import loadData
import allUserDefinedFuncs
import timeit
import pandas as pd


def haveDataM():
    return loadData.load_dataM()


def haveDataC():
    return loadData.load_dataC()


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
    refTaxa=None,
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

    allFunc = allUserDefinedFuncs.allUserFunc()
    results = []
    
    start = timeit.default_timer()
    stop = timeit.default_timer()
    print('Time: ', stop - start) 
    
CovData = haveDataC()
MicrobData = haveDataM()

IFAA(haveDataM(), haveDataC(), linkIDname = "id")
