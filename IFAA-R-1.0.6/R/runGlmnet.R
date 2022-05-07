#-------------------------------------------------------------------
## function for cross validation using glmnet package
#-------------------------------------------------------------------

runGlmnet=function(
  x,
  y,
  nPredics,
  family="gaussian",
  nfolds=10,
  lambda.min.ratio=0.05,
  nLam=100,
  intercept=TRUE,
  zeroSDCut=10^(-20)
){
  print("RunGlmnet Phase I")
  if (! "np" %in% ls() )   np <- reticulate::import("numpy")
  results=list()

  nBeta=ncol(x)

  # remove near constant x columns
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd,drop=FALSE]
  }
  rm(sdX)
  
  np$random$seed(as.integer(1))
  foldid=np$random$choice(as.integer(10), as.integer(length(y)), replace=T)+1
  foldid <- as.integer(foldid)
  
    cvResul=cv.glmnet(x=x,y=as.vector(y),nlambda=nLam,foldid = foldid,
                      family=family,intercept=intercept,standardize=FALSE)

    lamOpi=as.numeric(cvResul$lambda.min)

    rm(cvResul)

  finalLassoRun=glmnet(x=x,y=as.vector(y),lambda=lamOpi,family=family,
                       intercept=intercept,standardize=FALSE)
  rm(x,y)

  # convert to regular matrix format
  finalLassoRunBeta=as.vector(finalLassoRun$beta)

  # convert back to the full beta if there near constant x columns
  if(length(xWithNearZeroSd)>0){
    betaTrans=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                  unSelectList=sort(xWithNearZeroSd),newBetaNoInt=finalLassoRunBeta)
    beta=betaTrans$finalBeta
    rm(betaTrans)

    rm(xWithNearZeroSd)
    } else {
     beta=finalLassoRunBeta
     }
  rm(finalLassoRun)

  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  rm(beta)
  return(results)
}

