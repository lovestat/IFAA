# IFAA(dataM, dataC, "id", sequentialRun=TRUE)

results <- IFAA(MicrobData = dataM,
                CovData = dataC,
                linkIDname = "id",
                testCov = c("v1"),
                ctrlCov = c("v2", "v3"),
                paraJobs = 4,
                nRef = 60,
                seed = 1,
                sequentialRun=F)



results$sig_results
results$full_results
results$analysisResults$finalizedBootRefTaxon
results$analysisResults$allRefTaxNam
results$covriateNames
