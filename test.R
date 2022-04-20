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
                bootB = 30,
                sequentialRun=F)



results$sig_results
results$full_results
results$analysisResults$finalizedBootRefTaxon
results$analysisResults$allRefTaxNam
results$covriateNames


library(glmnet)
x = read.csv("../IFAA/x.csv")[, -1]
y = read.csv("../IFAA/y.csv")[, -1]

fit <- glmnet(as.matrix(x), y, 
              alpha=1, 
              standardize = FALSE, 
              intercept=TRUE)
fit$lambda
coef(fit, s = 0.01) %>% head(10)

cvfit <- glmnet::cv.glmnet(as.matrix(x), y, 
                           intercept = T, 
                           standardize = F)
cvfit$lambda
cvfit$lambda.min
coef(cvfit, s = cvfit$lambda.min)%>% head(10)







