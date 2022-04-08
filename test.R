# IFAA(dataM, dataC, "id", sequentialRun=TRUE)
results <- IFAA(MicrobData = dataM,
                CovData = dataC,
                linkIDname = "id",
                refTaxa = refTaxa[1:40],
                testCov = c("v1"),
                ctrlCov = c("v2","v3"),
                fdrRate = 0.15,
                # paraJobs = 4,
                seed = 1,
                sequentialRun=T)


library(tidyverse)

view(
  as.data.frame(as.matrix(xTildalong))
  
)

aaa = map(omegaRoot,  ~ as.data.frame(as.matrix(.x)))

bbb = aaa[[40]]

saveRDS(aaa, "aaa.rds")

results$sig_results
results$full_results



results$analysisResults$allRefTaxNam
results$covriateNames
names(results)
# [1] "covariatesData"  "covriateNames"   "analysisResults"
# [4] "sig_results"     "full_results"    "testCov"        
# [7] "ctrlCov"         "microbName"      "bootB"          
# [10] "refReadsThresh"  "balanceCut"      "SDThresh"       
# [13] "SDquantilThresh" "nRef"            "seed"           
# [16] "totalTimeMins"  

refTaxa = c("rawCount12", "rawCount31", "rawCount1", "rawCount2", "rawCount3", 
            "rawCount5", "rawCount6", "rawCount11", "rawCount13", "rawCount14", 
            "rawCount16", "rawCount20", "rawCount21", "rawCount22", "rawCount24", 
            "rawCount25", "rawCount26", "rawCount27", "rawCount29", "rawCount42", 
            "rawCount52", "rawCount33", "rawCount34", "rawCount44", "rawCount45", 
            "rawCount46", "rawCount54", "rawCount47", "rawCount48", "rawCount37", 
            "rawCount49", "rawCount50", "rawCount39", "rawCount51", "rawCount55", 
            "rawCount56", "rawCount53", "rawCount58", "rawCount59", "rawCount4", 
            "rawCount7", "rawCount10", "rawCount15", "rawCount17", "rawCount19", 
            "rawCount23", "rawCount28", "rawCount32", "rawCount35", "rawCount38", 
            "rawCount57", "rawCount40", "rawCount30", "rawCount8", "rawCount43", 
            "rawCount9", "rawCount60", "rawCount18", "rawCount41", "rawCount36"
)
