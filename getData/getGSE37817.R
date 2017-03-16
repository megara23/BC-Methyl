setwd("C:/Users/Windows User/Desktop/MethylMix/")
library(GEOquery)
dd = getGEO("GSE37817")
sapply(dd, annotation)

GSE37817.p.expr = pData(dd[[1]])
GSE37817.expr = exprs(dd[[1]])
GSE37817.p.methyl = pData(dd[[2]])
GSE37817.methyl = exprs(dd[[2]])

samples.match = colnames(GSE37817.expr) == rownames(GSE37817.p.expr)
if (!samples.match) {
  stop("need to match samples")
}

samples.match = colnames(GSE37817.methyl) == rownames(GSE37817.p.methyl)
if (!samples.match) {
  stop("need to match samples")
}

tmp = rep(NA, length(GSE37817.p.expr$title))
tmp[  grep("Control", GSE37817.p.expr$title)  ] = 0
tmp[  grep("Primary", GSE37817.p.expr$title)  ] = 1
GSE37817.expr.tumor = tmp

tmp2 = rep(NA, length(GSE37817.p.methyl$title))
tmp2[  grep("Control", GSE37817.p.methyl$title)  ] = 0
tmp2[  grep("Bladder cancer", GSE37817.p.methyl$title)  ] = 1
GSE37817.methyl.tumor = tmp2

save(GSE37817.p.methyl, GSE37817.methyl, GSE37817.methyl.tumor, file = "GSE37817.RData")


