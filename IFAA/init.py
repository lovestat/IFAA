#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 21:36:24 2022

@author: jin
"""


#!/usr/bin/python
# -*- coding: utf-8 -*-
from loadData import *
import allUserDefinedFuncs
from metaData import *
import timeit
import numpy as np
import pandas as pd
import warnings

maxDimensionScr=434*5*10**5
CovData = load_dataC()
MicrobData = load_dataM()
linkIDname = "id"
testCov = ["v1", "v2"]
ctrlCov = ["v3"]
testMany=True 
ctrlMany=False 
nRef=40 
nRefMaxForEsti=2 
refTaxa=[] 
adjust_method='BY' 
fdrRate=0.15 
paraJobs=None 
bootB=500 
standardize=False 
sequentialRun=False 
refReadsThresh=0.2 
taxkeepThresh=1 
SDThresh=0.05 
SDquantilThresh=0 
balanceCut=0.2 
seed=1 

i = 1