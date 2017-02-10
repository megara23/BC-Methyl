setwd("C:/Users/Windows User/Desktop/MethylMix/")
library(GEOquery)
dir = getGEO("GSE28094")

GSE28094.meth.p = pData(dir[[1]])
GSE28094.meth = exprs(dir[[1]])

GSE28094.meth.p = GSE28094.meth.p[GSE28094.meth.p$source_name_ch1 == 'bladder tumor tissue' |  GSE28094.meth.p$source_name_ch1 =='bladder normal tissue',] 
x = match(rownames(GSE28094.meth.p), colnames(GSE28094.meth))
GSE28094.meth = GSE28094.meth[,x]

tmp = rep(NA, length(GSE28094.meth.p$characteristics_ch1.4))
tmp[  grep("biomaterial: bladder normal tissue", GSE28094.meth.p$characteristics_ch1.4)  ] = 0
tmp[  grep("biomaterial: bladder tumor tissue", GSE28094.meth.p$characteristics_ch1.4)  ] = 1
GSE28094.meth.tumor = tmp


save(GSE28094.meth, GSE28094.meth.p, GSE28094.methyl.tumor, file = "GSE28094.RData")


        

