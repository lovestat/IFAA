#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 14:49:50 2022

@author: jin
"""

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