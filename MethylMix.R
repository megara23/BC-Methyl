load("C:/Users/Windows User/Desktop/MethylMix/GSE37817.RData")
library(MethylMix)
library(GEOquery)
MetNormal = GSE37817.methyl[,1:6]
MetCancer = GSE37817.methyl[,7:24]
ExpCancer = GSE37817.expr[,1:6]
data(MAcancer)
data(METcancer)
data(METnormal)

GPL = getGEO("GPL6102")
GPLmeth = getGEO("GPL8490")

genes.db = Table(GPL)
m = match(rownames(ExpCancer), genes.db$ID)
genes = genes.db$'ILMN_Gene'[m]
expresslist = data.frame(genes)
row.names(ExpCancer) = expresslist[1:length(rownames(ExpCancer)),]
ExpCancer = as.matrix(ExpCancer, na.rm=TRUE)


methyl.db = Table(GPLmeth)
m2 = match(rownames(MetCancer), methyl.db$ID)
cancermet = methyl.db$'Symbol' [m2]
cmethyl_list = data.frame(cancermet)
row.names(MetCancer) = cmethyl_list[1:length(rownames(MetCancer)),]
MetCancer = as.matrix(MetCancer, na.rm=TRUE)

m3 = match(rownames(MetNormal), methyl.db$ID)
normalmet = methyl.db$'Symbol' [m3]
nmethyl_list = data.frame(normalmet)
row.names(MetNormal) = nmethyl_list[1:length(rownames(MetNormal)),]
MetNormal = as.matrix(MetNormal, na.rm=TRUE)

MethylMix(MetCancer, MetNormal, ExpCancer)
