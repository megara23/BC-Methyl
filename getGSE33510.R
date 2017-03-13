setwd("C:/Users/Windows User/Desktop/MethylMix/")
library(GEOquery)
dir = getGEO("GSE33510")

sapply(dir, annotation)

GSE33510.meth.p = pData(dir[[1]])
GSE33510.meth.p = GSE33510.meth.p[GSE33510.meth.p$source_name_ch1 != 'Normal Ureter Sample',]
GSE33510.meth = exprs(dir[[1]])
x = match(rownames(GSE33510.meth.p), colnames(GSE33510.meth))
GSE33510.meth = GSE33510.meth[,x]
samples.match = colnames(GSE33510.meth) == rownames(GSE33510.meth.p)
if (!samples.match) {
  stop("need to match samples")
}

tmp = rep(NA, length(GSE33510.meth.p$title))
tmp[  grep("N", GSE33510.meth.p$title)  ] = 0
tmp[  grep("UC", GSE33510.meth.p$title)  ] = 1
GSE33510.meth.tumor_normal = tmp

save(GSE33510.meth, GSE33510.meth.p, GSE33510.meth.tumor_normal, file = "GSE33510.RData")
