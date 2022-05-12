# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:09:17 2022

@author: ss
"""

from IFAA import *
import numpy as np
from loadData import *


res = IFAA(load_dataM(), 
           load_dataC(), 
           testCov = ['v1'],
           ctrlCov = ['v2', 'v3'],
           paraJobs = 4,
           linkIDname="id",
           refTaxa = ["rawCount" + str(i + 1) for i in range(40)],
           bootB = 100,
           sequentialRun=True)

res['sig_results']
res['full_results']
res['analysisResults']['finalizedBootRefTaxon']
res['covriateNames']



## Test Binary Var
import numpy as np
from loadData import *

dataC = load_dataC()
dataC['v4'] = [5]*20 + [6]*20
dataC['v5'] = [5]*10 + [6]*30
dataC['v6'] = [0]*40
dataC['v7'] = [0]*40

res_bin = IFAA(load_dataM(), 
           dataC, 
           testCov = ['v1', 'v4', 'v5', 'v6', 'v7'],
           ctrlCov = ['v2', 'v3'],
           paraJobs = 6,
           linkIDname="id",
           refTaxa = ["rawCount" + str(i + 1) for i in range(40)],
           bootB = 100,
           sequentialRun=False)

res_bin['sig_results']
res_bin['full_results']
res_bin['analysisResults']['finalizedBootRefTaxon']
res_bin['covriateNames']



y = np.arange(1, 11)
y = np.array([7, 3, 10, 1, 2, 6, 8, 5, 4, 9])
x = np.arange(1, 1001).reshape((10, -1), order = 'F')
lm_res = OLS(y, x).fit()
lm_res.summary()






















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

len(glmnetCoef(fit, s = scipy.float64([0.01]), exact = False))

cvfit = cvglmnet(x = x.copy(), 
                 y = y.copy(), 
                 ptype = 'mse', 
                 nfolds = 10,
                 alpha = 1, 
                 standardize = False,
                 intr = True)
cvfit['lambda_min']

cvglmnetCoef(cvfit, s = 'lambda_min')

len(cvglmnetCoef(cvfit, s = 'lambda_min'))



import time

def timer(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        func(*args, **kwargs)

        print('The function ran for', time.time() - start)
    return wrapper

mat = np.random.normal(size=[20,18])

@timer
def np_svd(x):
    [np.linalg.svd(x) for _ in range(1000)]
    
@timer
def np_qr(x):
    [np.linalg.svd(x) for _ in range(1000)]

np_svd(mat)
np_qr(mat)


mat = np.array(
    [ [1,1,2],
      [1,2,4],
      [1,3,6],
      [1,4,8]]
    )

detectLDcol(mat)


mat = np.array(
    [ [1,1,1,2],
      [1,1,2,4],
      [1,1,3,6],
      [1,1,4,8]]
    )

detectLDcol(mat)

n = 50
inputs = np.random.normal(size = n*3).reshape(-1, 3)
x = np.column_stack( (inputs[:,2],
                      -.25 * inputs[:,2], 
                      inputs[:,0] + inputs[:,1],
                      inputs[:,0], 
                      inputs[:,1]) )
detectLDcol(x)

mat = x
