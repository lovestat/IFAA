#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 19:37:49 2022

@author: jin
"""

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
    results={}
    
    # get the original sample size
    nSub=nrow(data)
    MVarNamLength=len(Mprefix)
    
    # get taxa variable names
    microPositions = data.columns.str.startswith(Mprefix)
    nTaxa = len(which(microPositions))
    taxaNames = data.columns[microPositions]
    rm(microPositions)
    
    ## to add if qualifyreftax 
    if qualifyRefTax:
        qualifyData = data.loc[rowSums(data.loc[:,taxaNames]>0)>=2, :]
        w=qualifyData.loc[:,taxaNames]
        nSubQualif=nrow(qualifyData)
        taxaOverThresh=taxaNames[ colSums(w>0)>=(nSubQualif*refReadsThresh) ]
        if(len(taxaOverThresh)==0):
            print("There are no taxa with presence over the threshold:",refReadsThresh,". Try lower the reference taxa reads threshold.","\n")
    
        # check the sd threshold
        sdTaxaOverThresh=np.zeros(len(taxaOverThresh))
        
        for i in range(len(taxaOverThresh)):
            taxa_i=w.loc[:, taxaOverThresh[i]].to_numpy()
            if np.sum(taxa_i>0)>1:
                sdTaxaOverThresh[i]=np.std(taxa_i[ (taxa_i>0)], ddof=1)
        
        results['sdTaxa']=sdTaxaOverThresh

        TaxaOverSdThresh=taxaOverThresh[(sdTaxaOverThresh>=SDThresh)]
        if(len(TaxaOverSdThresh)==0):
            print("There are no taxa with SD over the SD threshold:",SDThresh, ". Try lower the SD threshold","\n")
            rm(taxa_i,taxaOverThresh)
    
        # check the sd quantile threshold
    
        sdAllTaxa=np.zeros(nTaxa)
        for i in range(nTaxa):
            taxaAll_i=w.loc[:,taxaNames[i]]
            posTaxaAll_i=taxaAll_i[(taxaAll_i>0)]
            if(len(posTaxaAll_i)>1):
                sdAllTaxa[i]=np.std(posTaxaAll_i, ddof=1)
                goodRefTaxaCandi=TaxaOverSdThresh[(sdAllTaxa>=np.quantile(sdAllTaxa,SDquantilThresh))]
                rm(sdAllTaxa,posTaxaAll_i,TaxaOverSdThresh)

        if len(goodRefTaxaCandi)==0:
            print("There are no taxa with SD over the SD quantile threshold:",SDquantilThresh, ". Try lower the SD quantile threshold","\n")
            rm(w)
    
    # get predictor data
    predNames=data.columns[data.columns.str.startswith(covsPrefix)].to_numpy()
    nPredics=len(predNames)
    
    ## to add if qualifyreftax 
    if qualifyRefTax:
        # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
        if len(binPredInd)>0:
            allBinPred=predNames[binPredInd]
            nBinPred=len(allBinPred)

            taxaBalanceBin=c()
            bin_nonz_sum = colSums(qualifyData.loc[:,allBinPred])
            
            if min(bin_nonz_sum,nSubQualif-bin_nonz_sum)<=np.floor(balanceCut*nSubQualif):
                raise Exception("one of the binary variable is not diverse enough")
            
            ## to add binary loop
            taxaBalanceBin = np.unique(taxaBalanceBin)
            # keep balanced taxa
            goodRefTaxaCandi=goodRefTaxaCandi[r_in(goodRefTaxaCandi, taxaBalanceBin)]
        results['goodRefTaxaCandi']=goodRefTaxaCandi
        rm(goodRefTaxaCandi)
    
    # return
    results['taxaNames']=taxaNames
    rm(taxaNames)
    results['predNames']=predNames
    rm(predNames)
    results['nTaxa']=nTaxa
    results['nSub']=nSub
    results['nPredics']=nPredics
    return results
    
    
    
    
