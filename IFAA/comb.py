#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 21:26:38 2022

@author: Shangchen
"""
import joblib
import contextlib
import math
import timeit
import scipy
import warnings
import numpy as np
import pandas as pd
import multiprocessing as mp
from joblib import Parallel, delayed
from tqdm import tqdm
from functools import partial
from sklearn.linear_model import Lasso
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
from statsmodels.api import OLS
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from itertools import compress
# import rpy2.robjects as robjects



def IFAA(
    MicrobData,
    CovData,
    linkIDname,
    testCov=np.empty(0),
    ctrlCov=np.empty(0),
    testMany=True,
    ctrlMany=False,
    nRef=40,
    nRefMaxForEsti=2,
    refTaxa=np.empty(0),
    adjust_method="fdr_by",
    fdrRate=0.15,
    paraJobs=np.empty(0),
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
    # Make arguments as numpy arries
    testCov = np.array(testCov)
    ctrlCov = np.array(ctrlCov)
    linkIDname = np.array(linkIDname)
    refTaxa = np.array(refTaxa)
    paraJobs = np.array([paraJobs])
    
    # results container, a dictionary
    results = {}

    start_time = timeit.default_timer()

    runMeta = metaData(
        MicrobData=MicrobData,
        CovData=CovData,
        linkIDname=linkIDname,
        taxkeepThresh=taxkeepThresh,
        testCov=testCov,
        ctrlCov=ctrlCov,
        testMany=testMany,
        ctrlMany=ctrlMany,
    )

    data = runMeta["data"]
    results["covariatesData"] = runMeta["covariatesData"]
    binaryInd = runMeta["binaryInd"]
    covsPrefix = runMeta["covsPrefix"]
    Mprefix = runMeta["Mprefix"]
    testCovInd = runMeta["testCovInd"]
    testCovInOrder = runMeta["testCovInOrder"]
    testCovInNewNam = runMeta["testCovInNewNam"]
    ctrlCov = runMeta["ctrlCov"]
    microbName = runMeta["microbName"]
    newMicrobNames = runMeta["newMicrobNames"]
    results["covriateNames"] = runMeta["xNames"]
    del runMeta

    if (refTaxa is not None) and (len(refTaxa) > 0):
        if sum(r_in(refTaxa, microbName)) != len(refTaxa):
            raise Exception(
                """
                             Error: One or more of the specified reference taxa in phase 1 have no sequencing reads 
                             or are not in the data set. Double check the names of the reference taxa and their 
                             sparsity levels."""
            )
    if nRefMaxForEsti < 2:
        nRefMaxForEsti = 2
        warnings.warn(
            "Warning: Needs at least two final reference taxon for estimation."
        )

    if nRef > len(microbName):
        raise Exception(
            "Error: number of random reference taxa can not be larger than the total number of taxa in the data. Try lower nRef"
        )

    refTaxa_newNam = newMicrobNames[r_in(microbName, refTaxa)]

    results["analysisResults"] = Regulariz(
        data=data,
        testCovInd=testCovInd,
        testCovInOrder=testCovInOrder,
        testCovInNewNam=testCovInNewNam,
        microbName=microbName,
        nRef=nRef,
        nRefMaxForEsti=nRefMaxForEsti,
        binaryInd=binaryInd,
        covsPrefix=covsPrefix,
        Mprefix=Mprefix,
        refTaxa=refTaxa_newNam,
        paraJobs=paraJobs,
        adjust_method=adjust_method,
        fwerRate=fdrRate,
        bootB=bootB,
        standardize=standardize,
        sequentialRun=sequentialRun,
        refReadsThresh=refReadsThresh,
        SDThresh=SDThresh,
        SDquantilThresh=SDquantilThresh,
        balanceCut=balanceCut,
        seed=seed,
    )
    
    rm(data)

    results['sig_results']=results['analysisResults']['sig_results']
    results['full_results']=results['analysisResults']['full_results']
    
    results['testCov']=testCovInOrder
    results['ctrlCov']=ctrlCov
    results['microbName']=microbName
    results['bootB']=bootB
    results['refReadsThresh']=refReadsThresh
    results['balanceCut']=balanceCut
    results['SDThresh']=SDThresh
    results['paraJobs']=paraJobs
    results['SDquantilThresh']=SDquantilThresh
    results['nRef']=nRef
    
    if isinstance(seed, int):
            results['seed']=seed
    else:
        results['seed']="No seed used."
                
                
    totalTimeMins =  (timeit.default_timer() - start_time)/60
    print("The entire analysis took ", np.round(totalTimeMins,2), " minutes")
                
    results['totalTimeMins']=totalTimeMins
                
    return results
            

def metaData(
    MicrobData,
    CovData,
    linkIDname,
    taxkeepThresh,
    testCov=np.empty(0),
    ctrlCov=np.empty(0),
    testMany=True,
    ctrlMany=True,
    MZILN=False,
):
    results = {}

    if not linkIDname:
        raise Exception("linkIDname is missing.")

    if len(testCov) > 0 and len(ctrlCov) > 0:
        if sum(r_in(np.concatenate((testCov, ctrlCov)), colnames(CovData))) != len(
            np.concatenate((testCov, ctrlCov))
        ):
            raise Exception("Error: some covariates are not available in the data.")

    if sum(r_in(testCov, ctrlCov)) > 0:
        warnings.warn(
            "Variables appeared in both testCov list and ctrlCov list will be treated as testCov."
        )

    # read microbiome data
    MdataWithId = MicrobData
    if len(colnames(MdataWithId)) != ncol(MdataWithId):
        raise Exception("Microbiome data lack variable names.")

    if (MdataWithId.loc[:, MdataWithId.columns != linkIDname].values < 0).any():
        raise Exception("Microbiome data contains negative values.")

    if MdataWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn(
            "There are over 80% missing values for the linkId variable in the Microbiome data file. Double check the data format."
        )

    # read covariate data
    CovarWithId = CovData
    if CovarWithId[linkIDname].isna().mean() > 0.8:
        warnings.warn(
            "There are over 80% missing values for the linkId variable in the covariates data file. Double check the data format."
        )

    Covariates1 = CovarWithId.loc[:, CovarWithId.columns != linkIDname]

    # determine testCov and ctrlCov
    if len(testCov) == 0:
        if not testMany:
            raise Exception(
                "No covariates are specified for estimating associations of interest."
            )
        else:
            print(
                "Associations are being estimated for all covariates since no covariates are specified for testCov."
            )
            testCov = colnames(Covariates1)
    results["testCov"] = testCov

    xNames = colnames(Covariates1)
    rm(Covariates1)

    if len(ctrlCov) == 0 and ctrlMany:
        print(
            "No control covariates are specified, all variables except testCov are considered as control covariates."
        )
        ctrlCov = xNames[r_ni(xNames, testCov)]
    ctrlCov = ctrlCov[r_ni(ctrlCov, testCov)]
    results["ctrlCov"] = ctrlCov

    # merge data to remove missing
    CovarWithId1 = CovarWithId[np.hstack([linkIDname, testCov, ctrlCov])]
    allRawData = CovarWithId1.merge(MdataWithId, on=linkIDname.tolist()).dropna()
    CovarWithId = allRawData.loc[:, r_in(colnames(allRawData), colnames(CovarWithId1))]
    Covariates = CovarWithId.loc[:, r_ni(colnames(CovarWithId), [linkIDname])]
    rm(CovarWithId1)

    if not all(Covariates.apply(lambda a: a.dtype.kind in "if")):
        raise Exception(
            "There are non-numeric variables in the covariates for association test."
        )

    MdataWithId = allRawData.loc[:, r_in(colnames(allRawData), colnames(MdataWithId))]
    Mdata_raw = MdataWithId.loc[:, r_ni(colnames(MdataWithId), linkIDname)]
    rm(allRawData)

    # check zero taxa and subjects with zero taxa reads
    numTaxaNoReads = sum(colSums(Mdata_raw != 0) <= taxkeepThresh)
    if numTaxaNoReads > 0:
        Mdata_raw = Mdata_raw.loc[:, colSums(Mdata_raw != 0) > taxkeepThresh]
        print(
            "There are ",
            numTaxaNoReads,
            " taxa without any sequencing reads before data merging, and excluded from the analysis",
        )
    rm(numTaxaNoReads)

    numSubNoReads = sum(rowSums(Mdata_raw != 0) <= 1)
    if numSubNoReads > 0:
        print(
            "There are ",
            numSubNoReads,
            " subjects with zero or one sequencing read and excluded from the analysis",
        )
        subKeep = inv_bool(rowSums(Mdata_raw != 0) <= 1)
        Mdata_raw = Mdata_raw.loc[subKeep, :]
        MdataWithId = MdataWithId.loc[subKeep, :]
        rm(subKeep)
    rm(numSubNoReads)

    Mdata = Mdata_raw
    rm(Mdata_raw)

    microbName1 = colnames(Mdata)
    microbName = microbName1
    newMicrobNames1 = np.array(["microb" + str(i + 1) for i in range(len(microbName))])
    newMicrobNames = newMicrobNames1
    results["Mprefix"] = "microb"
    Mdata = Mdata.rename(columns=dict(zip(microbName1, newMicrobNames)))

    MdataWithId_new = cbind([MdataWithId.loc[:, linkIDname], Mdata])

    results["microbName"] = microbName
    results["newMicrobNames"] = newMicrobNames

    if Covariates.isna().sum().sum() > 0:
        print("Samples with missing covariate values are removed from the analysis.")

    ## to add non-numeric add

    xNames = colnames(Covariates)
    nCov = len(xNames)

    ## to add binary check

    binCheck = Covariates.nunique()

    if any(binCheck == 2):
        binaryInd = which(binCheck == 2)
        results["varNamForBin"] = xNames[binCheck == 2]
        results["nBinVars"] = len(results["varNamForBin"])

        for i in range(results["varNamForBin"]):
            mini = min(Covariates.iloc[:, i])
            maxi = max(Covariates.iloc[:, i])
            if any(mini != 0 and maxi != 1):
                Covariates.iloc[
                    Covariates.iloc[:, i] == mini, results["varNamForBin"][i]
                ] = 0
                Covariates.iloc[
                    Covariates.iloc[:, i] == maxi, results["varNamForBin"][i]
                ] = 1
                print(
                    "Binary covariate",
                    i,
                    "is not coded as 0/1 which may generate analysis bias. It has been changed to 0/1. The changed covariates data can be extracted from the result file.",
                )

    results["nBinVars"] = 0
    binaryInd = []
    results["varNamForBin"] = []

    results["binaryInd"] = []
    results["xNames"] = colnames(Covariates)
    xNewNames = np.array(["x" + str(i + 1) for i in range(len(xNames))])
    Covariates = Covariates.rename(columns=dict(zip(colnames(Covariates), xNewNames)))
    results["covsPrefix"] = "x"
    results["xNewNames"] = xNewNames

    results["testCovInd"] = which(r_in(results["xNames"], testCov))

    results["testCovInOrder"] = results["xNames"][results["testCovInd"]]
    results["testCovInNewNam"] = results["xNewNames"][results["testCovInd"]]
    del (xNames, xNewNames)

    CovarWithId_new = cbind([CovarWithId.loc[:, linkIDname], Covariates])
    data = MdataWithId_new.merge(CovarWithId_new, on=linkIDname.tolist())
    dataOmit = data.dropna()

    results["covariatesData"] = CovarWithId_new
    results["covariatesData"].rename(
        columns=dict(
            zip(results["covariatesData"], np.insert(results["xNames"], 0, linkIDname))
        )
    )
    del (MdataWithId_new, CovarWithId_new)

    Mdata_omit = dataOmit.loc[:, np.array(newMicrobNames1)]

    # check taxa with zero or 1 read again after all missing data removed
    # to add

    results["data"] = dataOmit

    return results


def dataInfo(
    data,
    Mprefix,
    covsPrefix,
    binPredInd,
    refReadsThresh=None,
    SDThresh=None,
    SDquantilThresh=None,
    balanceCut=None,
    qualifyRefTax=False,
):
    results = {}

    # get the original sample size
    nSub = nrow(data)
    MVarNamLength = len(Mprefix)

    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    nTaxa = len(which(microPositions))
    taxaNames = data.columns[microPositions]
    rm(microPositions)

    ## to add if qualifyreftax
    if qualifyRefTax:
        qualifyData = data.loc[rowSums(data.loc[:, taxaNames] > 0) >= 2, :]
        w = qualifyData.loc[:, taxaNames]
        nSubQualif = nrow(qualifyData)
        taxaOverThresh = taxaNames[colSums(w > 0) >= (nSubQualif * refReadsThresh)]
        if len(taxaOverThresh) == 0:
            print(
                "There are no taxa with presence over the threshold:",
                refReadsThresh,
                ". Try lower the reference taxa reads threshold.",
                "\n",
            )

        # check the sd threshold
        sdTaxaOverThresh = np.zeros(len(taxaOverThresh))

        for i in range(len(taxaOverThresh)):
            taxa_i = w.loc[:, taxaOverThresh[i]].to_numpy()
            if np.sum(taxa_i > 0) > 1:
                sdTaxaOverThresh[i] = np.std(taxa_i[(taxa_i > 0)], ddof=1)

        results["sdTaxa"] = sdTaxaOverThresh

        TaxaOverSdThresh = taxaOverThresh[(sdTaxaOverThresh >= SDThresh)]
        if len(TaxaOverSdThresh) == 0:
            print(
                "There are no taxa with SD over the SD threshold:",
                SDThresh,
                ". Try lower the SD threshold",
                "\n",
            )
            rm(taxa_i, taxaOverThresh)

        # check the sd quantile threshold

        sdAllTaxa = np.zeros(nTaxa)
        for i in range(nTaxa):
            taxaAll_i = w.loc[:, taxaNames[i]]
            posTaxaAll_i = taxaAll_i[(taxaAll_i > 0)]
            if len(posTaxaAll_i) > 1:
                sdAllTaxa[i] = np.std(posTaxaAll_i, ddof=1)
                goodRefTaxaCandi = TaxaOverSdThresh[
                    (sdAllTaxa >= np.quantile(sdAllTaxa, SDquantilThresh))
                ]
                rm(sdAllTaxa, posTaxaAll_i, TaxaOverSdThresh)

        if len(goodRefTaxaCandi) == 0:
            print(
                "There are no taxa with SD over the SD quantile threshold:",
                SDquantilThresh,
                ". Try lower the SD quantile threshold",
                "\n",
            )
            rm(w)

    # get predictor data
    predNames = data.columns[data.columns.str.startswith(covsPrefix)].to_numpy()
    nPredics = len(predNames)

    ## to add if qualifyreftax
    if qualifyRefTax:
        # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
        if len(binPredInd) > 0:
            allBinPred = predNames[binPredInd]
            nBinPred = len(allBinPred)

            taxaBalanceBin = np.array([])
            bin_nonz_sum = colSums(qualifyData.loc[:, allBinPred])

            if min(bin_nonz_sum, nSubQualif - bin_nonz_sum) <= np.floor(
                balanceCut * nSubQualif
            ):
                raise Exception("one of the binary variable is not diverse enough")

            ## to add binary loop
            taxaBalanceBin = np.unique(taxaBalanceBin)
            # keep balanced taxa
            goodRefTaxaCandi = goodRefTaxaCandi[r_in(goodRefTaxaCandi, taxaBalanceBin)]
        results["goodRefTaxaCandi"] = goodRefTaxaCandi
        rm(goodRefTaxaCandi)

    # return
    results["taxaNames"] = taxaNames
    rm(taxaNames)
    results["predNames"] = predNames
    rm(predNames)
    results["nTaxa"] = nTaxa
    results["nSub"] = nSub
    results["nPredics"] = nPredics
    return results


def dataRecovTrans(data, ref, Mprefix, covsPrefix, xOnly=False, yOnly=False):
    results = {}

    # load A and log-ratio transformed RA
    data_and_init = AIcalcu(data=data, ref=ref, Mprefix=Mprefix, covsPrefix=covsPrefix)
    rm(data)

    taxaNames = data_and_init["taxaNames"]
    A = data_and_init["A"]
    logRatiow = data_and_init["logRatiow"]
    nSub = data_and_init["nSub"]
    nTaxa = data_and_init["nTaxa"]
    xData = data_and_init["xData"]
    nPredics = data_and_init["nPredics"]
    twoList = data_and_init["twoList"]
    lLast = data_and_init["lLast"]
    L = data_and_init["l"]
    lengthTwoList = data_and_init["lengthTwoList"]
    rm(data_and_init)

    nNorm = nTaxa - 1
    xDimension = nPredics + 1  # predictors+intercept

    # create omegaRoot
    omegaRoot = {}
    for j in range(lengthTwoList):
        i = twoList[j]
        if lLast[i] == nTaxa-1:
            omegaRoot[i] = np.eye(int(L[i] - 1))
        else:
            if L[i] == 2:
                omegaRoot[i] = np.sqrt(0.5)
            else:
                dim = L[i] - 1
                a = (1 + (dim - 2) / 2) / (1 / 2 * (1 + (dim - 1) / 2))
                b = -1 / (1 + (dim - 1) / 2)

                # calculate the square root of omega assuming it is exchangeable
                aStar = dim**2 / ((dim - 1) ** 2)
                bStar = b * (dim - 2) / (dim - 1) - a * (dim**2 - 2 * dim + 2) / (
                    (dim - 1) ** 2
                )
                cStar = (0.5 * b - 0.5 * a * (dim - 2) / (dim - 1)) ** 2
                cSquare = (-bStar + np.sqrt(bStar**2 - 4 * aStar * cStar)) / (
                    2 * aStar
                )

                if cSquare < 0:
                    raise Exception("no solution for square root of omega")
                    
                d = np.sqrt((0.5 * a - cSquare) / (dim - 1))

                if d is None:
                    raise Exception(
                        "no solution for off-diagnal elements for square root of omega"
                    )
                c = (0.5 * b - (dim - 2) * (d**2)) / (2 * d)
                omegaRoot[i] = -(
                    (c - d) * np.eye(int(dim)) + d * np.ones((int(dim), int(dim)))
                )

        rm(L, lLast)
        

    if xOnly:
        # create X_i in the regression equation using Kronecker product
        xDataWithInter = xData.copy()
        xDataWithInter.insert(0, "Intercept", 1)
        xDataWithInter = xDataWithInter.to_numpy()
        rm(xData)

        for j in range(lengthTwoList):
            i = twoList[j]
            xInRegres_i = np.kron(np.eye(nNorm), xDataWithInter[i, :])
            xDataTilda_i = omegaRoot[i] @ A[i] @ xInRegres_i
            rm(xInRegres_i)

            if j == 0:
                xTildalong = xDataTilda_i
            else:
                xTildalong = np.vstack((xTildalong, xDataTilda_i))

        rm(xDataWithInter, xDataTilda_i, omegaRoot, logRatiow)

        results["xTildalong"] = xTildalong
        rm(xTildalong)
        return results

    if yOnly:
        for j in range(lengthTwoList):
            i = twoList[j]
            Utilda_i = omegaRoot[i] @ logRatiow[i]
            if j == 0:
                UtildaLong = Utilda_i
            else:
                UtildaLong = np.hstack((UtildaLong, Utilda_i))

        rm(omegaRoot, logRatiow, Utilda_i)
        results["UtildaLong"] = UtildaLong
        rm(UtildaLong)
        return results


    # if not xOnly and yOnly
    # create X_i in the regression equation using Kronecker product
    xDataWithInter = xData.copy()
    rm(xData)
    xDataWithInter.insert(0, "Inter", 1)
    xDataWithInter = xDataWithInter.to_numpy()

    for j in range(lengthTwoList):
        i = twoList[j]
        xInRegres_i = np.kron(np.eye(nNorm), xDataWithInter[i, :])
        xDataTilda_i = omegaRoot[i] @ A[i] @ xInRegres_i
        rm(xInRegres_i)

        if j == 0:
            xTildalong = xDataTilda_i
        else:
            xTildalong = np.vstack((xTildalong, xDataTilda_i))

    rm(xDataWithInter, xDataTilda_i)

    for j in range(lengthTwoList):
        i = twoList[j]
        Utilda_i = omegaRoot[i] @ logRatiow[i]
        if j == 0:
            UtildaLong = Utilda_i
        else:
            UtildaLong = np.hstack((UtildaLong, Utilda_i))

    rm(omegaRoot, logRatiow, Utilda_i)

    # return objects
    results["UtildaLong"] = UtildaLong
    rm(UtildaLong)
    results["xTildalong"] = xTildalong
    rm(xTildalong)
    results["taxaNames"] = taxaNames
    rm(taxaNames)
    return results


def AIcalcu(data, ref, Mprefix, covsPrefix):
    results = {}
    # get the original sample size
    nSub = nrow(data)
    MVarNamLength = len(Mprefix)

    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    nTaxa = len(which(microPositions))
    nNorm = nTaxa - 1
    taxaNames = data.columns[microPositions]
    rm(microPositions)

    # rearrange taxa names
    otherTaxaNames = taxaNames[r_ni(taxaNames, ref)]
    taxaNames = np.hstack([otherTaxaNames, ref])

    # get predictor data
    predNames = data.columns[data.columns.str.startswith(covsPrefix)].to_numpy()
    nPredics = len(predNames)

    # taxa data
    w = data.loc[:, taxaNames]

    # extract x data
    xData = data.loc[:, predNames]
    rm(data, predNames)

    # transform data using log-ratio, creat Ai and Li
    l = np.empty(nSub)
    lLast = np.empty(nSub)
    taxa_non0 = {}
    taxa_0 = {}
    logRatiow = {}
    A = {}
    

    
    for i in range(nSub):
        taxa_nonzero = which(w.iloc[i, :] != 0)
        lLast[i] = np.max(taxa_nonzero)
        taxa_zero = which(w.iloc[i, :] == 0)
        taxa_non0[i] = w.iloc[i, taxa_nonzero]
        taxa_0[i] = w.iloc[i, taxa_zero]
        if len(taxa_nonzero) > 0:
            last_nonzero = np.max(taxa_nonzero)
            logwi = np.log(w.iloc[i, taxa_nonzero])
            l[i] = len(logwi)
            if l[i] > 1:
                logRatiow[i] = logwi[:-1:] - logwi[-1]
                zero_m = np.zeros((int(l[i]) - 1, nNorm))
                if last_nonzero == nTaxa-1:
                    aRow = np.arange(int(l[i]) - 1)
                    aCol = taxa_nonzero[:-1:]
                    zero_m[aRow, aCol] = 1
                else:
                    aRow = np.arange(int(l[i]) - 1)
                    aCol = taxa_nonzero[:-1:]
                    zero_m[aRow, aCol] = 1
                    zero_m[:, int(taxa_nonzero[int(l[i]) - 1])] = -1
                A[i] = zero_m
                rm(zero_m)
            else:
                logRatiow[i] = None
                A[i] = None
        else:
            l[i] = 0
            logRatiow[i] = None
            A[i] = None
            
    # obtain the list of samples whose have at least 2 non-zero taxa
    twoList = which(l > 1)
    lengthTwoList = len(twoList)

    rm(w)

    results["xData"] = xData
    rm(xData)

    results["logRatiow"] = logRatiow
    rm(logRatiow)
    results["A"] = A
    rm(A)
    results["twoList"] = twoList
    rm(twoList)
    results["taxaNames"] = taxaNames
    rm(taxaNames)
    results["lengthTwoList"] = lengthTwoList
    results["lLast"] = lLast
    results["l"] = l
    results["nTaxa"] = nTaxa
    results["nNorm"] = nNorm
    results["nSub"] = nSub
    results["nPredics"] = nPredics
    return results


def runScrParal(
    data,
    testCovInd,
    testCovInOrder,
    testCovInNewNam,
    nRef,
    paraJobs,
    refTaxa,
    standardize,
    sequentialRun,
    refReadsThresh,
    SDThresh,
    SDquantilThresh,
    balanceCut,
    Mprefix,
    covsPrefix,
    binPredInd,
    adjust_method,
    seed,
    maxDimensionScr=0.8 * 434 * 10 * 10**4,
):

    results = {}

    # load data info
    basicInfo = dataInfo(
        data=data,
        Mprefix=Mprefix,
        covsPrefix=covsPrefix,
        binPredInd=binPredInd,
        refReadsThresh=refReadsThresh,
        SDThresh=SDThresh,
        SDquantilThresh=SDquantilThresh,
        balanceCut=balanceCut,
        qualifyRefTax=True,
    )
    
    taxaNames = basicInfo["taxaNames"]
    nTaxa = basicInfo["nTaxa"]
    nPredics = basicInfo["nPredics"]
    nSub = basicInfo["nSub"]
    predNames = basicInfo["predNames"]

    results["goodRefTaxaCandi"] = basicInfo["goodRefTaxaCandi"]
    rm(basicInfo)

    nNorm = nTaxa - 1
    nAlphaNoInt = nPredics * nNorm
    nAlphaSelec = nPredics * nTaxa
    

    # make reference taxa list
    if len(refTaxa) < nRef:
        if seed is not None:
            np.random.seed(seed)

        taxon_to_be_sample = results["goodRefTaxaCandi"][
            r_ni(results["goodRefTaxaCandi"], refTaxa)
        ]
        num_to_be_sample = nRef - len(refTaxa)

        if num_to_be_sample >= len(taxon_to_be_sample):
            num_to_be_sample = len(taxon_to_be_sample)
        print(
            "The number of candidate reference taxon is smaller than the number of taxon required in phase 1. The number of taxon was set to be ",
            num_to_be_sample,
        )

        refTaxa_extra = np.random.choice(
            taxon_to_be_sample, num_to_be_sample, replace=False
        )
        refTaxa = np.hstack((refTaxa, refTaxa_extra))

        if len(refTaxa) == 0:
            raise Exception(
                "No candidate reference taxon is available. Please try to lower the reference taxon boundary."
            )

    if len(refTaxa) >nRef:
        if seed is not None:
            np.random.seed(seed)
        refTaxa = np.random.choice(refTaxa, nRef, replace=True)
    
    results["refTaxa"] = np.array(refTaxa)

    ## run original data screen
    screen1 = originDataScreen(
        data=data,
        testCovInd=testCovInd,
        nRef=nRef,
        refTaxa=refTaxa,
        paraJobs=paraJobs,
        Mprefix=Mprefix,
        covsPrefix=covsPrefix,
        binPredInd=binPredInd,
        standardize=standardize,
        sequentialRun=sequentialRun,
        adjust_method=adjust_method,
        seed=seed,
    )

    results["countOfSelecForAPred"] = screen1["countOfSelecForAPred"]
    results["estOfSelectForAPred"] = screen1["estOfSelectForAPred"]
    results["testCovCountMat"] = screen1["testCovCountMat"]
    results["testEstMat"] = screen1["testEstMat"]

    # rm(screen1)

    nTestCov = len(testCovInd)
    results["nTestCov"] = nTestCov
    results["nTaxa"] = nTaxa
    results["nPredics"] = nPredics

    results["taxaNames"] = taxaNames
    # rm(taxaNames)
    return results


def dataSparsCheck(data, Mprefix):
    results = {}

    # get the original sample size
    nSub = nrow(data)
    MVarNamLength = len(Mprefix)

    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    taxaNames = data.columns[microPositions]
    rm(microPositions)

    w = data.loc[:, taxaNames]
    rm(data, taxaNames)
    overallSparsity = np.round(100 * np.mean(w.values == 0), 2)
    print(overallSparsity, "percent of microbiome sequencing reads are zero")

    # check zero taxa and subjects with zero taxa reads
    numTaxaNoReads = sum(colSums(w) == 0)
    if numTaxaNoReads > 0:
        print(
            "There are ",
            numTaxaNoReads,
            " taxa without any sequencing reads and excluded from the analysis",
        )
    rm(numTaxaNoReads)

    numSubNoReads = sum(rowSums(w) == 0)
    if numSubNoReads > 0:
        print(
            "There are ",
            numSubNoReads,
            " subjects without any sequencing reads and excluded from the analysis.",
        )
    rm(numSubNoReads, w)


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
    refReadsThresh,
    SDThresh,
    SDquantilThresh,
    balanceCut,
    adjust_method,
    seed,
):

    results = {}
    regul_start_time = timeit.default_timer()

    nTestCov = len(testCovInd)
    dataSparsCheck(data=data, Mprefix=Mprefix)

    # load abundance data info

    binCheck = data.loc[:, testCovInNewNam].nunique()
    binaryInd = which(binCheck == 2)

    data.info = dataInfo(
        data=data, Mprefix=Mprefix, covsPrefix=covsPrefix, binPredInd=binaryInd
    )
    nSub = data.info["nSub"]
    taxaNames = data.info["taxaNames"]
    nPredics = data.info["nPredics"]
    nTaxa = data.info["nTaxa"]
    rm(data.info)

    regul_start_time = timeit.default_timer()
    print("Start Phase 1 analysis")

    selectRegroup = getScrResu(
        data=data,
        testCovInd=testCovInd,
        testCovInOrder=testCovInOrder,
        testCovInNewNam=testCovInNewNam,
        nRef=nRef,
        paraJobs=paraJobs,
        refTaxa=refTaxa,
        standardize=standardize,
        sequentialRun=sequentialRun,
        refReadsThresh=refReadsThresh,
        SDThresh=SDThresh,
        SDquantilThresh=SDquantilThresh,
        balanceCut=balanceCut,
        Mprefix=Mprefix,
        covsPrefix=covsPrefix,
        binPredInd=binaryInd,
        adjust_method=adjust_method,
        seed=seed,
    )
    nRef_smaller = np.max((2, math.ceil(nRef / 2)))
    while_loop_ind = False
    loop_num = 0
    print("33 percent of phase 1 analysis has been done")
    while while_loop_ind is False:
        loop_num = loop_num + 1
        if loop_num >= 2:
            break
        if len(refTaxa) < nRef_smaller:
            refTaxa_smaller = np.hstack(
                (
                    refTaxa,
                    (selectRegroup["goodIndpRefTaxWithCount"].index)[
                        : nRef_smaller - len(refTaxa)
                    ],
                )
            )
        else:
            refTaxa_smaller = refTaxa
        fin_ref_1 = selectRegroup["finalIndpRefTax"]
        ref_taxa_1 = selectRegroup["refTaxa"]
        selectRegroup = getScrResu(
            data=data,
            testCovInd=testCovInd,
            testCovInOrder=testCovInOrder,
            testCovInNewNam=testCovInNewNam,
            nRef=nRef_smaller,
            paraJobs=paraJobs,
            refTaxa=refTaxa_smaller,
            standardize=standardize,
            sequentialRun=sequentialRun,
            refReadsThresh=refReadsThresh,
            SDThresh=SDThresh,
            SDquantilThresh=SDquantilThresh,
            balanceCut=balanceCut,
            Mprefix=Mprefix,
            covsPrefix=covsPrefix,
            binPredInd=binaryInd,
            adjust_method=adjust_method,
            seed=seed,
        )
        fin_ref_2 = selectRegroup["finalIndpRefTax"]
        ref_taxa_2 = selectRegroup["refTaxa"]
        while_loop_ind = np.array_equal(fin_ref_1, fin_ref_2) | np.array_equal(
            ref_taxa_1, ref_taxa_2
        )

        if not while_loop_ind:
            print(
                np.round(100 * (loop_num + 1) / 3, 0),
                " percent of phase 1 analysis has been done",
            )
        if while_loop_ind:
            print("100 percent of phase 1 analysis has been done")



    results["selecCountOverall"] = selectRegroup["selecCountOverall"]
    results["selecCountOverall"].columns = microbName
    results["selecCountMatIndv"] = selectRegroup["selecCountMatIndv"]
    finalIndpRefTax = microbName[r_in(taxaNames, selectRegroup["finalIndpRefTax"])]
    results["finalRefTaxonQualified"] = selectRegroup["refTaxonQualified"]
    results["goodIndpRefTaxLeastCount"] = microbName[
        r_in(taxaNames, selectRegroup["goodIndpRefTaxLeastCount"])
    ]
    results["goodIndpRefTaxWithCount"] = selectRegroup["goodIndpRefTaxWithCount"]
    results["goodIndpRefTaxWithCount"].index = microbName[
        np.hstack(
            [
                which(r_in(taxaNames, i))
                for i in selectRegroup["goodIndpRefTaxWithCount"].index
            ]
        )
    ]

    results["goodIndpRefTaxWithEst"] = selectRegroup["goodIndpRefTaxWithEst"]
    results["goodIndpRefTaxWithEst"].index = microbName[
        np.hstack(
            [
                which(r_in(taxaNames, i))
                for i in selectRegroup["goodIndpRefTaxWithEst"].index
            ]
        )
    ]

    results["goodRefTaxaCandi"] = microbName[
        r_in(taxaNames, selectRegroup["goodRefTaxaCandi"])
    ]
    results["randomRefTaxa"] = microbName[r_in(taxaNames, selectRegroup["refTaxa"])]
    goodIndpRefTax_ascend = results["goodIndpRefTaxWithCount"].sort_values()
    goodIndpRefTaxNam = goodIndpRefTax_ascend.index
    rm(selectRegroup)

    MCPExecuTime = (timeit.default_timer() - regul_start_time) / 60
    results["MCPExecuTime"] = MCPExecuTime
    print("Phase 1 analysis used ", np.round(MCPExecuTime, 2), " minutes")

    results["finalizedBootRefTaxon"] = finalIndpRefTax

    startT = timeit.default_timer()
    print("Start Phase 2 parameter estimation")

    unestimableTaxa = np.empty(0)
    qualifyData = data

    binCheck = qualifyData.nunique()
    binaryInd = which(binCheck == 2)

    # to add
    if len(binaryInd) > 0:
        pass

    # check zero taxa and subjects with zero taxa reads
    TaxaNoReads = which(colSums(qualifyData.loc[:, taxaNames]) == 0).astype(int)
    rm(qualifyData)
    unestimableTaxa = np.unique(np.hstack((unestimableTaxa, taxaNames[TaxaNoReads])))
    results["unEstTaxa"] = microbName[r_in(taxaNames, unestimableTaxa)]
    
    ## numpy unique will sort it !! avoid it by pandas
    allRefTaxNam = pd.unique(
        np.concatenate((results["finalizedBootRefTaxon"], goodIndpRefTaxNam))
    )
    nGoodIndpRef = len(allRefTaxNam)
    results["allRefTaxNam"] = allRefTaxNam

    results["nRefUsedForEsti"] = np.min((nGoodIndpRef, nRefMaxForEsti))

    results["estiList"] = {}

    for iii in range(results["nRefUsedForEsti"]):
        print("Start estimation for the ", iii, "th final reference taxon")
        time11 = timeit.default_timer()
        originTaxNam = allRefTaxNam[iii]
        newRefTaxNam = taxaNames[r_in(microbName, originTaxNam)]
        results["estiList"][originTaxNam] = bootResuHDCI(
            data=data,
            refTaxa=newRefTaxNam,
            originRefTaxNam=originTaxNam,
            bootB=bootB,
            binPredInd=binaryInd,
            covsPrefix=covsPrefix,
            Mprefix=Mprefix,
            testCovInOrder=testCovInOrder,
            adjust_method=adjust_method,
            microbName=microbName,
            fwerRate=fwerRate,
            paraJobs=paraJobs,
            standardize=standardize,
            seed=seed,
        )
        time12 = timeit.default_timer()
        print(
            "Estimation done for the ",
            iii,
            "th final reference taxon and it took ",
            round((time12 - time11) / 60, 3),
            " minutes",
        )

    endT = timeit.default_timer()

    print(
        "Phase 2 parameter estimation done and took ",
        round((endT - startT) / 60, 3),
        " minutes.",
    )
    
    
    ### calculate mean ###
    fin_ref_taxon_name=np.array(list(results['estiList'].keys()))
    
    all_cov_sig_list={}
    all_cov_list={}
    
    for i in range(len(fin_ref_taxon_name)):
        i_name= fin_ref_taxon_name[i]
        all_cov_list[i_name]=results['estiList'][i_name]['all_cov_list']
        all_cov_sig_list[i_name]=results['estiList'][i_name]['sig_list_each']
    
    results['all_cov_list']=all_cov_list
    results['all_cov_sig_list']=all_cov_sig_list

    
    
    ref_taxon_name=np.array(list( all_cov_list.keys() ))
    
    all_cov_list_0 = all_cov_list[ref_taxon_name[0]]
    all_cov_list_1 = all_cov_list[ref_taxon_name[1]]
        
    exclu_0=r_ni( colnames(all_cov_list_0['est_save_mat']), ref_taxon_name[1] )
    exclu_1=r_ni( colnames(all_cov_list_1['est_save_mat']), ref_taxon_name[0] )
    
    est_save_mat_mean=(all_cov_list_0['est_save_mat'].loc[:,exclu_0]+all_cov_list_1['est_save_mat'].loc[:,exclu_1])/2
    se_mat_mean=(all_cov_list_0['se_mat'].loc[:,exclu_0]+all_cov_list_1['se_mat'].loc[:,exclu_1])/2
    CI_low_mat_mean=(all_cov_list_0['CI_low_mat'].loc[:,exclu_0]+all_cov_list_1['CI_low_mat'].loc[:,exclu_1])/2
    CI_up_mat_mean=(all_cov_list_0['CI_up_mat'].loc[:,exclu_0]+all_cov_list_1['CI_up_mat'].loc[:,exclu_1])/2
    
    ## to add
    # if len(binaryInd)>0:
    #     est_save_mat_mean.loc[biNoneryInd,results.loc['unEstTaxa']]=None
    #     CI_low_mat_mean.loc[biNoneryInd,results.loc['unEstTaxa']]=None
    #     CI_up_mat_mean.loc[biNoneryInd,results.loc['unEstTaxa']]=None
    #     se_mat_mean.loc[biNoneryInd,results.loc['unEstTaxa']]=None
    
    p_value_unadj_mean = (est_save_mat_mean/se_mat_mean).apply(
        lambda x: (1-scipy.stats.norm.cdf(np.abs(x)))*2, 
        )
        
    p_value_adj_mean = p_value_unadj_mean.apply(
        lambda x: multipletests(x, method=adjust_method)[1],
        axis = 1,
        result_type='expand'
        )
    
    p_value_adj_mean.columns = p_value_unadj_mean.columns
    
    colname_use=colnames(est_save_mat_mean)
    
    
    sig_ind = np.where((p_value_adj_mean < fwerRate).to_numpy())
    sig_row = sig_ind[0]
    sig_col = sig_ind[1]
    
    est_sig=[est_save_mat_mean.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    CI_low_sig=[CI_low_mat_mean.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    CI_up_sig=[CI_up_mat_mean.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    p_adj_sig=[p_value_adj_mean.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    se_sig=[se_mat_mean.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    
    est_sig=np.array(est_sig)
    CI_low_sig=np.array(CI_low_sig)
    CI_up_sig=np.array(CI_up_sig)
    p_adj_sig=np.array(p_adj_sig)
    se_sig=np.array(se_sig)
    
    cov_sig_index=np.sort(np.unique(sig_row))

    sig_list_each_mean={}
    
    if len(cov_sig_index)>0:
        for iii in range(len(cov_sig_index)):
            sig_loc=which(sig_row==cov_sig_index[iii]).astype(int)
            est_spe_cov=est_sig[sig_loc]
            CI_low_spe_cov=CI_low_sig[sig_loc]
            CI_up_spe_cov=CI_up_sig[sig_loc]
            p_adj_spe_cov=p_adj_sig[sig_loc]
            se_spe_cov=se_sig[sig_loc]
            
            cov_sig_mat=np.zeros( (len(sig_loc),5) )
            cov_sig_mat = pd.DataFrame(cov_sig_mat, 
                                       columns = ["estimate","SE est","CI low","CI up","adj p-value"])
            cov_sig_mat.iloc[:,0]=est_spe_cov
            cov_sig_mat.iloc[:,1]=se_spe_cov
            cov_sig_mat.iloc[:,2]=CI_low_spe_cov
            cov_sig_mat.iloc[:,3]=CI_up_spe_cov
            cov_sig_mat.iloc[:,4]=p_adj_spe_cov
            
            cov_sig_mat.index=colname_use[sig_col[sig_loc]]
            sig_list_each_mean[ testCovInOrder[cov_sig_index[iii]] ] = cov_sig_mat
    
    results['sig_results']=sig_list_each_mean
    full_results={}
    
    for j in range(len(testCovInOrder)):
        j_name = testCovInOrder[j]
        est_res_save_all=pd.concat((est_save_mat_mean.loc[j_name,],
                               se_mat_mean.loc[j_name,],
                               CI_low_mat_mean.loc[j_name,],
                               CI_up_mat_mean.loc[j_name,],
                               p_value_adj_mean.loc[j_name,]),
                                   axis = 1)
        est_res_save_all.columns= ["estimate","SE est","CI low","CI up","adj p-value"]
        full_results[j_name]=est_res_save_all
        
    results['full_results']=full_results

    results['nTaxa']=nTaxa
    results['nPredics']=nPredics
    
    # return results

    results['nRef']=nRef
    return results
    

def bootResuHDCI(
    data,
    refTaxa,
    originRefTaxNam,
    bootB,
    binPredInd,
    covsPrefix,
    Mprefix,
    testCovInOrder,
    adjust_method,
    microbName,
    fwerRate,
    paraJobs,
    standardize,
    seed,
    maxDimension=434 * 5 * 10 ** 5,
    bootLassoAlpha=0.05,
):
    

    results = {}
    
    # load data info
    basicInfo = dataInfo(
        data=data, Mprefix=Mprefix, covsPrefix=covsPrefix, binPredInd=binPredInd
    )

    taxaNames = basicInfo["taxaNames"]
    ii=which( r_in(basicInfo['taxaNames'], refTaxa))
    nTaxa = basicInfo["nTaxa"]
    nPredics = basicInfo["nPredics"]
    rm(basicInfo)

    nNorm = nTaxa - 1
    nAlphaNoInt = nPredics * nNorm
    nAlphaSelec = nPredics * nTaxa

    countOfSelec = np.zeros(nAlphaSelec)
    
    resultsByRefTaxon={}
    
    # inital Lasso OLS estimate
    dataForEst=dataRecovTrans(data=data,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix)

    x=dataForEst['xTildalong']
    y=dataForEst['UtildaLong']
    rm(dataForEst)

    xCol=x.shape[1]
    
    maxSubSamplSiz=np.min( (50000,math.floor(maxDimension/xCol))  )
    nToSamplFrom=len(y)
    subSamplK=math.ceil(nToSamplFrom/maxSubSamplSiz)
    if subSamplK==1:
        maxSubSamplSiz=nToSamplFrom
    
    nRuns=(math.ceil(subSamplK/3))
    
    if x.shape[0] > x.shape[1] :
        for k in range(nRuns):
            rowToKeep = np.random.choice(nToSamplFrom, maxSubSamplSiz, replace=False)
            xSub = x[rowToKeep, :]
            ySub = y[rowToKeep]
            
            lm_res = OLS(ySub, xSub).fit()
            
            lm_res.summary()
            lm_res.params
            lm_res.bse
            lm_res.tvalues
            lm_res.pvalues
            
            bootResu = np.transpose(
                np.vstack( (lm_res.params,
                lm_res.bse,
                lm_res.tvalues,
                lm_res.pvalues) ))

            if k==0 :
                bootResu_k =  bootResu
            else:
                bootResu_k=bootResu_k+bootResu
        

        fin_ref_taxon_name=originRefTaxNam
        all_cov_list={}
        nTestcov=len(testCovInOrder)
        
        boot_est=bootResu_k[:,0]/nRuns
        se_est_all=bootResu_k[:,1]/nRuns
        boot_CI=lm_res.conf_int()
        
        ref_taxon_name=originRefTaxNam
        
        p_value_save_mat=np.zeros( (nTestcov, nTaxa-1) )
        est_save_mat=np.zeros( (nTestcov, nTaxa-1) )
        CI_up_mat=np.zeros( (nTestcov, nTaxa-1) )
        CI_low_mat=np.zeros( (nTestcov, nTaxa-1) )
        se_mat=np.zeros( (nTestcov, nTaxa-1) )
        

        
        for ii in range(nTestcov):
            se_est=se_est_all[(ii+1):len(lm_res.params):(nPredics+1)]
            boot_est_par=boot_est[(ii+1):len(lm_res.params):(nPredics+1)]
            
            p_value_unadj=(1-scipy.stats.norm.cdf(np.abs(boot_est_par/se_est)))*2
            boot_est_CI_low=boot_est_par-1.96*se_est
            boot_est_CI_up=boot_est_par+1.96*se_est
            
            p_value_adj=multipletests(
                p_value_unadj, method=adjust_method
            )[1]
            
            p_value_save_mat[ii,:]=p_value_adj
            est_save_mat[ii,]=boot_est_par
            CI_low_mat[ii,]=boot_est_CI_low
            CI_up_mat[ii,]=boot_est_CI_up
            se_mat[ii,]=se_est
            
    else:
        for k in range(nRuns):
            pass
            ## to add HDCI
            
    colname_use=microbName[microbName!=ref_taxon_name]
    
    p_value_save_mat=pd.DataFrame(p_value_save_mat, index = testCovInOrder, columns = colname_use)
    est_save_mat=pd.DataFrame(est_save_mat, index = testCovInOrder, columns = colname_use)
    CI_low_mat=pd.DataFrame(CI_low_mat, index = testCovInOrder, columns = colname_use)
    CI_up_mat=pd.DataFrame(CI_up_mat, index = testCovInOrder, columns = colname_use)
    se_mat=pd.DataFrame(se_mat, index = testCovInOrder, columns = colname_use)
    
    sig_ind = np.where((p_value_save_mat < fwerRate).to_numpy())
    sig_row = sig_ind[0]
    sig_col = sig_ind[1]
    
    est_sig=[est_save_mat.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    CI_low_sig=[CI_low_mat.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    CI_up_sig=[CI_up_mat.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    p_adj_sig=[p_value_save_mat.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    se_sig=[se_mat.iloc[sig_row[my_i], sig_col[my_i]] for my_i in range(len(sig_row))]
    
    est_sig=np.array(est_sig)
    CI_low_sig=np.array(CI_low_sig)
    CI_up_sig=np.array(CI_up_sig)
    p_adj_sig=np.array(p_adj_sig)
    se_sig=np.array(se_sig)
    
    cov_sig_index=np.sort(np.unique(sig_row))

    sig_list_each={}
    
    if len(cov_sig_index)>0:
        for iii in range(len(cov_sig_index)):
            sig_loc=which(sig_row==cov_sig_index[iii]).astype(int)
            est_spe_cov=est_sig[sig_loc]
            CI_low_spe_cov=CI_low_sig[sig_loc]
            CI_up_spe_cov=CI_up_sig[sig_loc]
            p_adj_spe_cov=p_adj_sig[sig_loc]
            se_spe_cov=se_sig[sig_loc]
            
            cov_sig_mat=np.zeros( (len(sig_loc),5) )
            cov_sig_mat = pd.DataFrame(cov_sig_mat, 
                                       columns = ["estimate","SE est","CI low","CI up","adj p-value"])
            cov_sig_mat.iloc[:,0]=est_spe_cov
            cov_sig_mat.iloc[:,1]=se_spe_cov
            cov_sig_mat.iloc[:,2]=CI_low_spe_cov
            cov_sig_mat.iloc[:,3]=CI_up_spe_cov
            cov_sig_mat.iloc[:,4]=p_adj_spe_cov
            
            cov_sig_mat.index=colname_use[sig_col[sig_loc]]
            sig_list_each[ testCovInOrder[cov_sig_index[iii]] ] = cov_sig_mat
            
            
    all_cov_list['est_save_mat']=est_save_mat
    all_cov_list['p_value_save_mat']=p_value_save_mat
    all_cov_list['CI_low_mat']=CI_low_mat
    all_cov_list['CI_up_mat']=CI_up_mat
    all_cov_list['se_mat']=se_mat
    
    results['sig_list_each']=sig_list_each
    results['all_cov_list']=all_cov_list
    
    return results

def getScrResu(
    data,
    testCovInd,
    testCovInOrder,
    testCovInNewNam,
    nRef,
    paraJobs,
    refTaxa,
    standardize,
    sequentialRun,
    refReadsThresh,
    SDThresh,
    SDquantilThresh,
    balanceCut,
    Mprefix,
    covsPrefix,
    binPredInd,
    adjust_method,
    seed,
    goodIndeCutPerc=0.33,
):
    results = {}
    # run permutation
    scrParal = runScrParal(
        data=data,
        testCovInd=testCovInd,
        testCovInOrder=testCovInOrder,
        testCovInNewNam=testCovInNewNam,
        nRef=nRef,
        paraJobs=paraJobs,
        refTaxa=refTaxa,
        standardize=standardize,
        sequentialRun=sequentialRun,
        refReadsThresh=refReadsThresh,
        SDThresh=SDThresh,
        SDquantilThresh=SDquantilThresh,
        balanceCut=balanceCut,
        Mprefix=Mprefix,
        covsPrefix=covsPrefix,
        binPredInd=binPredInd,
        adjust_method=adjust_method,
        seed=seed,
    )
    

    selecCountOverall = scrParal["countOfSelecForAPred"]
    selecEstOverall = scrParal["estOfSelectForAPred"]

    selecCountMatIndv = scrParal["testCovCountMat"]
    selecEstMatIndv = scrParal["testEstMat"]

    taxaNames = scrParal["taxaNames"]
    goodRefTaxaCandi = scrParal["goodRefTaxaCandi"]

    nTaxa = scrParal["nTaxa"]
    nPredics = scrParal["nPredics"]
    nTestCov = scrParal["nTestCov"]
    results["refTaxa"] = scrParal["refTaxa"]
    rm(scrParal)

    if nTestCov == 1:
        results["selecCountMatIndv"] = selecCountOverall
        results["selecEstMatIndv"] = selecEstOverall
    if nTestCov > 1:
        results["selecCountMatIndv"] = selecCountMatIndv
        results["selecEstMatIndv"] = selecEstMatIndv
        rm(selecCountMatIndv)

    goodIndpRefTaxWithCount = selecCountOverall.iloc[
        0, r_in(colnames(selecCountOverall), goodRefTaxaCandi)
    ]
    goodIndpRefTaxWithEst = selecEstOverall.iloc[
        0, r_in(colnames(selecEstOverall), goodRefTaxaCandi)
    ]

    if len(goodIndpRefTaxWithCount) == 0:
        results["goodIndpRefTaxLeastCount"] = np.array([])
    else:
        results["goodIndpRefTaxLeastCount"] = goodIndpRefTaxWithCount.index[
            np.lexsort((np.abs(goodIndpRefTaxWithEst), goodIndpRefTaxWithCount))
        ][0:2]
        goodIndpRefTaxWithEst = np.abs(
            goodIndpRefTaxWithEst[
                np.lexsort((np.abs(goodIndpRefTaxWithEst), goodIndpRefTaxWithCount))
            ]
        )
        goodIndpRefTaxWithCount = goodIndpRefTaxWithCount[
            np.lexsort((np.abs(goodIndpRefTaxWithEst), goodIndpRefTaxWithCount))
        ]

    results["selecCountOverall"] = selecCountOverall
    results["goodIndpRefTaxWithCount"] = goodIndpRefTaxWithCount
    results["goodIndpRefTaxWithEst"] = goodIndpRefTaxWithEst
    results["goodRefTaxaCandi"] = goodRefTaxaCandi
    rm(goodRefTaxaCandi)
    results["refTaxonQualified"] = 2
    results["finalIndpRefTax"] = results["goodIndpRefTaxLeastCount"]

    return results


def originDataScreen(
    data,
    testCovInd,
    nRef,
    paraJobs,
    refTaxa,
    standardize,
    sequentialRun,
    Mprefix,
    covsPrefix,
    binPredInd,
    adjust_method,
    seed,
    maxDimensionScr=434 * 5 * 10**5,
):

    results = {}

    # load data info
    basicInfo = dataInfo(
        data=data, Mprefix=Mprefix, covsPrefix=covsPrefix, binPredInd=binPredInd
    )

    taxaNames = basicInfo["taxaNames"]
    nTaxa = basicInfo["nTaxa"]
    nPredics = basicInfo["nPredics"]
    rm(basicInfo)

    nNorm = nTaxa - 1
    nAlphaNoInt = nPredics * nNorm
    nAlphaSelec = nPredics * nTaxa

    countOfSelec = np.zeros(nAlphaSelec)

    # overwrite nRef if the reference taxon is specified
    nRef = len(refTaxa)

    forEachUnitRun_partial = partial(
        forEachUnitRun,
        taxaNames,
        refTaxa,
        Mprefix,
        covsPrefix,
        maxDimensionScr,
        nPredics,
        data,
        nAlphaSelec,
        nAlphaNoInt,
        nTaxa,
        standardize
    )

    startT1 = timeit.default_timer()
    if len(paraJobs) == 0:
        availCores = mp.cpu_count()
        if isinstance(availCores, int):
            paraJobs = max(1, availCores - 2)

    if not sequentialRun:
        print(
            paraJobs,
            " parallel jobs are registered for analyzing ",
            nRef,
            " reference taxa in Phase 1",
        )

        with tqdm_joblib(tqdm(desc="Progress", total=nRef)) as progress_bar:
            scr1Resu = Parallel(n_jobs=int(paraJobs))(
                delayed(forEachUnitRun_partial)(i) for i in range(nRef)
            )

    if sequentialRun:
        print(
            " Sequential running analysis for ",
            nRef,
            " reference taxa in Phase 1",
        )
        
        scr1Resu = [forEachUnitRun_partial(i) for i in tqdm(range(nRef), desc="Sequential")]

    endT = timeit.default_timer()
    

    scr1ResuSelec = np.hstack([i["selection"][:, np.newaxis] for i in scr1Resu])
    scr1ResuEst = np.hstack([i["coef"][:, np.newaxis] for i in scr1Resu])

    # create count of selection for individual testCov
    countOfSelecForAllPred = scr1ResuSelec.sum(axis=1).reshape((-1, nPredics))
    countOfSelecForAllPred = np.transpose(countOfSelecForAllPred)
    EstOfAllPred = scr1ResuEst.sum(axis=1).reshape((-1, nPredics))
    EstOfAllPred = np.transpose(EstOfAllPred)

    testCovCountMat = countOfSelecForAllPred[
        testCovInd, :
    ]
    testEstMat = EstOfAllPred[
        testCovInd, :
    ]
    # rm(scr1ResuSelec, testCovInd, countOfSelecForAllPred, EstOfAllPred)
    
    
    # create overall count of selection for all testCov as a whole
    countOfSelecForAPred = testCovCountMat.sum(axis=0).reshape((1, -1))
    estOfSelectForAPred = testEstMat.sum(axis=0).reshape((1, -1))

    countOfSelecForAPred = pd.DataFrame(countOfSelecForAPred)
    estOfSelectForAPred = pd.DataFrame(estOfSelectForAPred)

    countOfSelecForAPred.columns = taxaNames
    estOfSelectForAPred.columns = taxaNames

    # return results
    results["testCovCountMat"] = testCovCountMat
    results["testEstMat"] = testEstMat
    # rm(testCovCountMat, testEstMat)
    results["countOfSelecForAPred"] = countOfSelecForAPred
    results["estOfSelectForAPred"] = estOfSelectForAPred
    # rm(countOfSelecForAPred, estOfSelectForAPred)
    return results


def forEachUnitRun(
    taxaNames,
    refTaxa,
    Mprefix,
    covsPrefix,
    maxDimensionScr,
    nPredics,
    data,
    nAlphaSelec,
    nAlphaNoInt,
    nTaxa,
    standardize,
    i,
):

    ii = which(taxaNames == refTaxa[i])
    dataForEst = dataRecovTrans(
        data=data, ref=refTaxa[i], Mprefix=Mprefix, covsPrefix=covsPrefix
    )

    xTildLongTild_i = dataForEst["xTildalong"]
    yTildLongTild_i = dataForEst["UtildaLong"]
    rm(dataForEst)

    maxSubSamplSiz = np.min(
        (50000.0, np.floor(maxDimensionScr / xTildLongTild_i.shape[1]))
    ).astype(int)

    nToSamplFrom = xTildLongTild_i.shape[0]

    subSamplK = np.ceil(nToSamplFrom / maxSubSamplSiz).astype(int)

    if subSamplK == 1:
        maxSubSamplSiz = nToSamplFrom

    nRuns = np.ceil(subSamplK / 3).astype(int)

    for k in range(nRuns):
        rowToKeep = np.random.choice(nToSamplFrom, maxSubSamplSiz, replace=False)

        x = xTildLongTild_i[rowToKeep, :]
        y = yTildLongTild_i[rowToKeep]
        
        
        if x.shape[0] > (3 * x.shape[1]):
            Penal_i = runlinear(x=x, y=y, nPredics=nPredics)
            BetaNoInt_k = (Penal_i["betaNoInt"] != 0).astype(int)
            EstNoInt_k = np.abs(Penal_i["coef_est_noint"])
        else:
            Penal_i = runGlmnet(x=x, y=y, nPredics=nPredics, standardize=standardize)
            BetaNoInt_k = (Penal_i["betaNoInt"] != 0).astype(int)
            EstNoInt_k = np.abs(Penal_i["betaNoInt"])

        rm(Penal_i)
        if k == 0:
            BetaNoInt_i = BetaNoInt_k
            EstNoInt_i = EstNoInt_k
        if k > 0:
            BetaNoInt_i = BetaNoInt_i + BetaNoInt_k
            EstNoInt_i = EstNoInt_i + EstNoInt_k
        rm(BetaNoInt_k, EstNoInt_k)

    rm(k, x, y, xTildLongTild_i)
    BetaNoInt_i = BetaNoInt_i / nRuns
    EstNoInt_i = EstNoInt_i / nRuns
    selection_i = np.zeros(nAlphaSelec)
    coef_i = np.zeros(nAlphaSelec)
        
        
    if ii == 0:
        np_assign_but(selection_i, np.linspace(0, nPredics - 1, nPredics), BetaNoInt_i)
        np_assign_but(coef_i, np.linspace(0, nPredics - 1, nPredics), EstNoInt_i)
    if ii == (nTaxa-1) :
        np_assign_but(
            selection_i,
            np.linspace(nAlphaSelec - nPredics, nAlphaSelec-1, nPredics),
            BetaNoInt_i,
        )
        np_assign_but(
            coef_i,
            np.linspace(nAlphaSelec - nPredics, nAlphaSelec-1, nPredics),
            EstNoInt_i,
        )
    if (ii > 0) & (ii < (nTaxa-1)):
        selection_i[0 : int((nPredics * (ii)))] = BetaNoInt_i[
            0 : int((nPredics * (ii)))
        ]
        selection_i[int(nPredics * (ii+1)) : (nAlphaSelec)] = BetaNoInt_i[
            int(nPredics * (ii)) : (nAlphaNoInt+1)
        ]
        coef_i[0 : int((nPredics * (ii)))] = EstNoInt_i[
            0 : int((nPredics * (ii)))
        ]
        coef_i[int(nPredics * (ii+1)) : (nAlphaSelec)] = EstNoInt_i[
            int(nPredics * (ii)) : (nAlphaNoInt+1)
        ]
    rm(BetaNoInt_i)
    # create return vector
    recturnlist = {}
    recturnlist["selection"] = selection_i
    recturnlist["coef"] = coef_i
    return recturnlist


def runlinear(x, y, nPredics, fwerRate=0.25, adjust_method="fdr_bh"):
    
    x = pd.DataFrame(x, columns = ["x" + str(i+1) for i in range(x.shape[1])] )
    y = pd.Series(y, name = "y")
    results = {}
    nBeta = x.shape[1]
    nObsAll = len(y)
    ## Build the data and formula
    dat = cbind( (y, x) )
    formula = dat.columns[0] + " ~ " + " + ".join(dat.columns[1:])  + " - 1"
    ## Fit the lm model
    lm_res = smf.ols(formula, data = dat).fit()
        
    p_value_est = lm_res.pvalues
    disc_index = np.arange(0, len(p_value_est), (nPredics + 1))
    p_value_est_noint = np.delete(p_value_est.to_numpy(), disc_index, axis=0)

    ## this method automatically convert over 1 values to 1
    p_value_est_noint_adj = multipletests(
        p_value_est_noint, alpha=0.05, method=adjust_method
    )[1]
    
    
    coef_est = np.abs(lm_res.params)
    disc_index = np.arange(0, len(p_value_est), (nPredics + 1))
    ## NA coef is snot considered here
    coef_est_noint = np.delete(coef_est.to_numpy(), disc_index, axis=0)

    # return
    results["betaNoInt"] = p_value_est_noint_adj < fwerRate
    # results["betaInt"] = p_value_est
    results["coef_est_noint"] = coef_est_noint

    return results


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
    zeroSDCut=10 ** (-20),
):
    results = {}
    nBeta = x.shape[1]
    nObsAll = len(y)

    # remove near constant x columns
    sdX = np.std(x, axis=0)
    xWithNearZeroSd = which(sdX <= zeroSDCut)
    if len(xWithNearZeroSd) > 0:
        x = np.delete(x, xWithNearZeroSd, axis=1)
    rm(sdX)

    # calculate lambda max
    if family == "gaussian":
        lamMax = np.max(np.abs((x * y[:, np.newaxis]).sum(0))) / nObsAll

    lamVec = np.linspace(lamMax, 0, nLam + 1)[0:nLam]
    if standardize is True:
        scaler = StandardScaler()
        x = scaler.fit_transform(x)

    cvStartTime = timeit.default_timer()
    cvStartTimei = timeit.default_timer()

    cvResul = LassoCV(
        alphas=lamVec, fit_intercept=intercept, cv=nfolds, n_jobs=-1, selection="random"
    ).fit(x, y)
    lamOpi = cvResul.alpha_

    cvExeTimei = (timeit.default_timer() - cvStartTimei) / 60
    cvAllTime = (timeit.default_timer() - cvStartTime) / 60

    finalLassoRun = Lasso(alpha=lamOpi, fit_intercept=intercept).fit(x, y)

    finalLassoRunBeta = finalLassoRun.coef_

    finalLassoRunBeta

    # convert back to the full beta if there near constant x columns
    if len(xWithNearZeroSd) > 0:
        pass
    else:
        beta = finalLassoRunBeta

    disc_index = np.arange(0, len(beta), (nPredics + 1))
    results["betaNoInt"] = np.delete(beta, disc_index, axis=0)

    return results


# =============================================================================
# Help Functions
# =============================================================================

def colnames(x):
    if not isinstance(x, pd.core.frame.DataFrame):
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.columns.to_numpy(dtype="U")


def ncol(x):
    if not isinstance(x, pd.core.frame.DataFrame):
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return len(x.columns)


def nrow(x):
    if not isinstance(x, pd.core.frame.DataFrame):
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return len(x.index)


def colSums(x):
    if not isinstance(x, pd.core.frame.DataFrame):
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis=0)


def rowSums(x):
    if not isinstance(x, pd.core.frame.DataFrame):
        raise Exception("Input is not pandas.core.frame.DataFrame ")
    return x.sum(axis=1)


def r_in(x, y):
    x, y = np.array(x), np.array(y)
    return np.array([np.isin(item, y) for item in x]).astype(bool)


def r_ni(x, y):
    x, y = np.array(x), np.array(y)
    return np.array([not np.isin(item, y) for item in x]).astype(bool)


def rm(*args):
    del args


def slice_bool(x, y):
    """Slicing list by a boolean list"""
    return list(compress(x, y))


def inv_bool(x):
    return [not i for i in x]


def cbind(x):
    return pd.concat(x, axis=1)


def which(x):
    return np.array([i for i, j in enumerate(x) if j])


def np_assign_but(ar, but_ind, v):
    assign_ind = np.setdiff1d(np.arange(len(ar)), but_ind)
    ar[assign_ind] = v
    return


@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""

    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()