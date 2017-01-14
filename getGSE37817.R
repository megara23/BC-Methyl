setwd("C:/Users/Windows User/Desktop/MethylMix/")
library(GEOquery)
dd = getGEO("GSE37817")
sapply(dd, annotation)

GSE37817.p.methyl = pData(dd[[2]])
GSE37817.methyl = exprs(dd[[2]])

samples.match = colnames(GSE37817.methyl) == rownames(GSE37817.p.methyl)
if (!samples.match) {
  stop("need to match samples")
}

tmp2 = rep(NA, length(GSE37817.p.methyl$title))
tmp2[  grep("Control", GSE37817.p.methyl$title)  ] = 0
tmp2[  grep("Bladder cancer", GSE37817.p.methyl$title)  ] = 1
GSE37817.methyl.tumor = tmp2

save(GSE37817.p.methyl, GSE37817.methyl, 
     GSE37817.methyl.tumor, file = "GSE37817.RData")
