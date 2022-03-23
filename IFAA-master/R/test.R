lapply(list.files("~/Dropbox (UFL)/Project-PhD/withDrLi/IFAA/IFAA-master/R"),
       function(x) source(glue::glue("~/Dropbox (UFL)/Project-PhD/withDrLi/IFAA/IFAA-master/R/{x}")))
load("~/Dropbox (UFL)/Project-PhD/withDrLi/IFAA/IFAA-master/data/dataC.rda")
load("~/Dropbox (UFL)/Project-PhD/withDrLi/IFAA/IFAA-master/data/dataM.rda")
pacman::p_load(methods,future,Matrix,HDCI,qlcMatrix,
               expm,rlecuyer,mathjaxr,glmnet,stats,utils,foreach,
               tidyverse)

MicrobData = dataM
CovData = dataC
linkIDname = "id"
taxkeepThresh = 0
testCov = c("v1", "v2")
ctrlCov = "v3"
ssfunc::eval_text(str_replace_all("ctrlCov=NULL,
  testMany=TRUE,
  ctrlMany=FALSE,
  nRef=40,
  nRefMaxForEsti=2,
  refTaxa=NULL,
  adjust_method='BY',
  fdrRate=0.15,
  paraJobs=NULL,
  bootB=500,
  standardize=FALSE,
  sequentialRun=FALSE,
  refReadsThresh=0.2,
  taxkeepThresh=1,
  SDThresh=0.05,
  SDquantilThresh=0,
  balanceCut=0.2,
  seed=1", ",", ";"))
