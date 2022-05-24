# ifaa
[![PyPI version](https://badge.fury.io/py/ifaa.svg)](https://badge.fury.io/py/ifaa)


## Overview

`ifaa` is a robust approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem. It can be also directly applied to relative abundance (RA) data to make inference on AA because the ratio of two RA is equal ratio of their AA. This algorithm can estimate and test the associations of interest while adjusting for potential confounders. High-dimensional covariates are handled with regularization. The estimates of this method have easy interpretation like a typical regression analysis. High-dimensional covariates are handled with regularization and it is implemented by parallel computing. False discovery rate is automatically controlled by this approach. Zeros do not need to be imputed by a positive value for the analysis. The IFAA package also offers the 'MZILN' function for estimating and testing associations of abundance ratios with covariates.

The following mixed effect model is used for the association: 

$\log(\mathcal{Y}_i^k)|\mathcal{Y}_i^k>0=\beta^{0k}+X_i^T\beta^k+W_i^T\gamma^k+Z_i^Tb_i+\epsilon_i^k,\hspace{0.2cm}k=1,...,K+1$

where

- $\mathcal{Y}_i^k$ is the AA of taxa $k$ in subject $i$ in the entire ecosystem. 

- $X_i$ is the covariate of interest

- $W_i$ is the confounder.

- $Z_i$ is the design matrix for random effects. 

- $\epsilon_i^k$ is the random error. 

- $\beta^k$ is the regression coefficients that will be estimated and tested with the `IFAA()` function.

The challenge in microbiome analysis is that we can not oberve $\mathcal{Y}_i^k$. What is observed is its small proportion: $Y_i^k=C_i\mathcal{Y}^k_i$ where $C_i$ is an unknown number between 0 and 1 that denote the observed proportion. The IFAA method successfuly addressed this challenge.

## Package installation 

To install, type the following command in Python console:


```python
pip install ifaa
```

The package could be also installed from GitHub using the following code: 


```python
pip install git+https://github.com/lovestat/ifaa
```

## Input for IFAA() function

Most of the time, users just need to feed the first five inputs to the function: `MicrobData`, `CovData`, `linkIDname`, `testCov` and `ctrlCov`. All other inputs can just take their default values. Below are all the inputs of the functions

- `MicrobData`: Microbiome data matrix containing microbiome absolute abundance or relative abundance with each row per sample and each column per taxon/OTU/ASV (or any other unit). The data matrix can have zero-valued data points. It should also contain an id variable to be linked with the id variable in the covariates data: `CovData`. This argument can also take file directory path. For example, `MicrobData="C://...//microbiomeData.tsv"`.

- `CovData`: Covariates data matrix containing covariates and confounders with each row per sample and each column per variable. Any categorical variable should be converted into dummy variables in this data matrix unless it can be treated as a continuous variable. It should also contain an id variable to be linked with the id variable in the microbiome data: `MicrobData`. This argument can also take file directory path. For example, `CovData="C://...//covariatesData.tsv"`.

- `linkIDname`: The common variable name of the id variable in both `MicrobData` and `CovData`. The two data sets will be merged by this id variable.

- `testCov`: Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.

- `ctrlCov`: Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.

- `testMany`: This takes logical value `True` or `False`. If `True`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `True` which does not do anything if `testCov` is not `NULL`.

- `ctrlMany`: This takes logical value `True` or `False`. If `True`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `False`.

- `nRef`: The number of randomly picked reference taxa used in phase 1. Default number is `40`. 
- `nRefMaxForEsti`: The maximum number of final reference taxa used in phase 2. The default is `2`.
- `refTaxa`: A vector of taxa names. These are reference taxa specified by the user to be used in phase 1 if the user believe these taxa are indepenent of the covariates. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly. 
- `adjust_method` The adjusting method used for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
- `fdrRate`: The false discovery rate for identifying taxa/OTU/ASV associated with `testCov`. Default is `0.15`.
- `paraJobs`: If `sequentialRun` is `False`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`. 
- `bootB`: Number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. The default is `500`.
- `standardize` This takes a logical value `True` or `False`. If `True`, the design matrix for X will be standardized in the analyses and the results. Default is `False`.
- `sequentialRun`: This takes a logical value `True` or `False`. Default is `False`. This argument could be useful for debug.
- `refReadsThresh`: The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
- `taxDropThresh`: The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
- `SDThresh`: The threshold of standard deviations of sequencing reads for been chosen as the reference taxon in phase 2. The default is `0.05` which means the standard deviation of sequencing reads should be at least `0.05` in order to be chosen as reference taxon.
- `SDquantilThresh`: The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as reference taxon. The default is `0`.
- `balanceCut`: The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default number is `0.2` which means at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
- `seed`: Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding. 

## Output for IFAA() function

The estimation results are saved in the following lists: 

- `sig_results`: A list containing estimating results that are statistically significant. 

- `full_results`: A list containing all estimating results. NA denotes unestimable.

The covariates data used in the analyses including `testCov` and `ctrlCov` is saved in the following object:

- `covariatesData`: A dataset containing covariates and confounders used in the analyses


## Examples

The example datasets `dataM` and `dataC` are included in the package. They could be accessed by: 

```python
import numpy as np
from loadData import *
```
Both the microbiome data `dataM` and the covariates data `dataC` contain 40 samples (i.e., 40 rows). 

- `dataM` contains 60 taxa with absolute abundances and these are gut microbiome.  

- `dataC` contains 3 covariates. 

Next we analyze the data to test the association between microbiome and the variable `"v1"` while adjusting for the variables (potential confounders) `"v2"` and `"v3"`.


```python
>>> res = IFAA(load_dataM().iloc[:,:],
...            load_dataC().iloc[:,:],
...            testCov = ['v1'],
...            ctrlCov = ['v2', 'v3'],
...            paraJobs = 4,
...            linkIDname="id",
...            refTaxa = ["rawCount" + str(i + 1) for i in range(40)],
...            bootB = 100,
...            sequentialRun=False)
```


In this example, we are only interested in testing the associations with `"v1"` which is why `testCov=c("v1")`. The variables `"v2" and "v3"` are adjusted as potential confounders in the analyses. The final analysis results are saved in the list `sig_results`: 

```python
res['sig_results']
# {'v1':             estimate    SE est    CI low     CI up   adj p-value
#  rawCount18  0.028144  0.005027  0.018292  0.037997  1.935964e-06
#  rawCount36  0.029991  0.005009  0.020174  0.039809  2.873697e-07
#  rawCount41  0.033031  0.005040  0.023152  0.042909  1.512894e-08}
```
The results found three taxa `"rawCount18"`, `"rawCount36"`, `"rawCount41"` associated with `"v1"` while adjusting for `"v2" and "v3"`. The regression coefficients and their 95% confidence intervals are provided. These coefficients correspond to $\beta^k$ in the model equation. 

The interpretation is that 

- Every unit increase in `"v1"` is associated with approximately 2.8% increase in the absolute abundance of `"rawCount18"`, approximately 2.9% increase in the absolute abundance of `"rawCount36"`, and approximately 3.3% increase in the absolute abundance of `"rawCount41"` in the entire gut ecosystem. 


All the analyzed covariates including `testCov` and `ctrlCov` can be extracted using the object `covariatesData`. The covariates data of the first 10 subjects can extracted as follows:


```python
res_bin['covariatesData'].iloc[0:10, :]
#    id          x1  x2  x3  x4  x5          x6          x7
# 0   1   58.069691   0   0   0   0  -49.903757  -15.306430
# 1   2   25.965216   0   0   0   0  -68.588941  -23.109922
# 2   3  193.716251   0   0   0   0  124.401856  119.567468
# 3   4   72.156467   0   0   0   0  -98.485358    2.877972
# 4   5   98.062712   0   0   0   0   23.553576  -79.893161
# 5   6   83.094848   0   0   0   0 -116.958209 -107.641285
# 6   7    8.217154   0   0   0   0 -205.644800 -139.958481
# 7   8   36.169820   0   0   0   0   58.957085   26.890379
# 8   9  152.786131   0   0   0   0  162.609349  138.731954
# 9  10   41.621790   0   0   0   0   65.154267   59.974310
```

## References 

- Zhigang Li, Lu Tian, A. James O'Malley, Margaret R. Karagas, Anne G. Hoen, Brock C. Christensen, Juliette C. Madan, Quran Wu, Raad Z. Gharaibeh, Christian Jobin, Hongzhe Li (2020) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. arXiv:1909.10101v3

- Zhigang Li, Katherine Lee, Margaret Karagas, Juliette Madan, Anne Hoen, James O’Malley and Hongzhe Li (2018 ) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences  10(3):587-608
