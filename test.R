# IFAA(dataM, dataC, "id", sequentialRun=TRUE)
library(tidyverse)

results <- IFAA(MicrobData = dataM[, ],
                CovData = dataC[, ],
                linkIDname = "id",
                testCov = c("v1"),
                ctrlCov = c("v2", "v3"),
                paraJobs = 4,
                refTaxa = glue::glue("rawCount{1:40}"),
                seed = 1,
                bootB = 100,
                sequentialRun=F)



results$sig_results
results$full_results
results$analysisResults$finalizedBootRefTaxon
results$analysisResults$allRefTaxNam
results$covriateNames



results.bin <- IFAA(MicrobData = dataM[, ],
                    CovData = cbind(dataC,
                                    data.frame(v4 = rep(5:6, each = 20),
                                               v5 = c(rep(5, 10), rep(6, 30) ))
                    ),
                    linkIDname = "id",
                    testCov = c("v1", "v4", "v5"),
                    ctrlCov = c("v2", "v3"),
                    paraJobs = 4,
                    refTaxa = glue::glue("rawCount{1:40}"),
                    seed = 1,
                    bootB = 100,
                    sequentialRun=F)



results.bin$sig_results
results.bin$full_results
results.bin$analysisResults$finalizedBootRefTaxon
results.bin$analysisResults$allRefTaxNam
results.bin$covriateNames


library(glmnet)
x = read.csv("../IFAA/x.csv")[, -1]
y = read.csv("../IFAA/y.csv")[, -1]

fit <- glmnet(as.matrix(x), y, 
              alpha=1, 
              standardize = FALSE, 
              intercept=TRUE)

fit$lambda
coef(fit, s = 0.0045) %>% head(10)
length(coef(fit, s = 0.0045))
cvfit <- glmnet::cv.glmnet(as.matrix(x), y, 
                           intercept = F, 
                           standardize = F)
cvfit$lambda
cvfit$lambda.min
coef(cvfit, s = cvfit$lambda.min)
length(coef(cvfit, s = cvfit$lambda.min))






