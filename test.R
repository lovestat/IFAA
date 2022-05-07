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

{
CovData = cbind(dataC,
                v4 = rep(5:6, each = 20),
                v5 = c(rep(5, 10), rep(6, 30)),
                v6 = 0,
                v7 = 0,
                v8 = rep(5:6, each = 20))

results.bin <- IFAA(MicrobData = dataM,
                    CovData = CovData,
                    linkIDname = "id",
                    testCov = c("v1", "v4", "v5", "v6", "v7"),
                    ctrlCov = c("v2", "v3"),
                    paraJobs = 4,
                    refTaxa = glue::glue("rawCount{1:40}"),
                    seed = 1,
                    bootB = 100,
                    sequentialRun=T)

}


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


## Zero SD variable
## v6 = all 0
## v7 = all 0


x = data.frame(x1 = 1:10,
               x2 = 2*(1:10),
               x3 = 3*(1:10),
               x4 = 1:10,
               x5 = rep(5, 10),
               x6 = rep(10, 10))

x$x4 = 2*x$x1
x$x5 = x$x1+x$x2
x$x6 = rnorm(10)

y = sample(10)
fit = lm(y~as.matrix(x)-1)
summary(fit)

fit$coefficients
summary(fit)$coefficients
