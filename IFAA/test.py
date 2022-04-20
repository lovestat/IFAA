# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:09:17 2022

@author: ss
"""

import numpy as np
from loadData import *

res = IFAA(load_dataM().iloc[:,:], 
           load_dataC().iloc[:,:], 
           testCov = ['v1'],
           ctrlCov = ['v2', 'v3'],
           paraJobs = 4,
           linkIDname="id",
           refTaxa = ["rawCount" + str(i + 1) for i in range(40)],
           bootB = 30,
           sequentialRun=True,
           CILearnerNam="Lasso")

res['sig_results']
res['full_results']
res['analysisResults']['finalizedBootRefTaxon']
res['covriateNames']

x = np.load("x.npy")
y = np.load("y.npy")

x_pd = pd.DataFrame(x)
y_pd = pd.DataFrame(y)

x_pd.to_csv("x.csv")
y_pd.to_csv("y.csv")

from sklearn import linear_model
from sklearn.model_selection import KFold


m = linear_model.Lasso(alpha = 0.01,
                       fit_intercept=True)
m.fit(x, y)
m.coef_[0:10]

kf = KFold(n_splits=10, shuffle = True)

cvm = linear_model.LassoCV(n_alphas=100,
                           fit_intercept=True,
                           cv = kf)
cvm.fit(x,y)
cvm.alpha_
cvm.coef_

cvm = linear_model.LassoLarsCV(fit_intercept=True,
                                cv = kf)

cvm.fit(x,y)
cvm.alpha_
cvm.coef_[0:10]

lm = linear_model.LinearRegression(fit_intercept=True)
lm.fit(x,y)
lm.coef_[0:10]

import glmnet_python
import scipy, importlib, pprint, matplotlib.pyplot as plt, warnings
from glmnet import glmnet; from glmnetPlot import glmnetPlot 
from glmnetPrint import glmnetPrint; from glmnetCoef import glmnetCoef; from glmnetPredict import glmnetPredict
from cvglmnet import cvglmnet; from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot; from cvglmnetPredict import cvglmnetPredict


fit = glmnet(x = x.copy(), 
             y = y.copy(), 
             family = 'gaussian', 
             alpha = 1, 
             standardize = False)

glmnetCoef(fit, s = scipy.float64([0.01]), exact = False)


cvfit = cvglmnet(x = x.copy(), 
                 y = y.copy(), 
                 ptype = 'mse', 
                 nfolds = 10,
                 alpha = 1, 
                 standardize = False)
cvfit['lambda_min']

cvglmnetCoef(cvfit, s = 'lambda_min')









