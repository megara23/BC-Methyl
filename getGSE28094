setwd("C:/Users/Windows User/Desktop/MethylMix/")
library(GEOquery)
dir = getGEO("GSE28094")

sapply(dir, annotation)

GSE28094.meth.p = pData(dir[[1]])
GSE28094.meth = exprs(dir[[1]])

samples.match = colnames(GSE28094.meth) == rownames(GSE28094.meth.p)
if (!samples.match) {
  stop("need to match samples")
}

keep1 = GSE28094.meth.p[,8] == "bladder tumor tissue"
keep2 = GSE28094.meth.p[,8] == "bladder normal tissue"
keep.all = keep1 | keep2

keep = function(x)
{
  x = readline(x)
  return(as.list(keep.all == TRUE))
}

wekeep = as.list(keep(GSE28094.meth.p))
wekeep = as.vector(wekeep == TRUE)
GSE28094.meth.p = GSE28094.meth.p[c(wekeep),]
GSE28094.meth = GSE28094.meth[,c(wekeep)]

tmp = rep(NA, nrow(GSE28094.meth.p))
tmp[GSE28094.meth.p[,8] == "bladder tumor tissue"] = 1
tmp[GSE28094.meth.p[,8] == "bladder normal tissue"]= 0
GSE28094.methyl.tumor = tmp

save(GSE28094.meth, GSE28094.meth.p, GSE28094.methyl.tumor, file = "GSE28094.RData")
