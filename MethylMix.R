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
reallykeep = expresslist[,1] != ""
row.names(ExpCancer) = expresslist[1:length(rownames(ExpCancer)),]
ExpCancer = as.matrix(ExpCancer, na.rm = TRUE)
ExpCancer = ExpCancer[reallykeep,]

methyl.db = Table(GPLmeth)
m2 = match(rownames(MetCancer), methyl.db$ID)
cancermet = methyl.db$'Symbol' [m2]
cmethyl_list = data.frame(cancermet)
keep = cmethyl_list[,1] != ""
row.names(MetCancer) = cmethyl_list[1:length(rownames(MetCancer)),]
MetCancer = as.matrix(MetCancer, na.rm = TRUE)
MetCancer = MetCancer[keep,]


m3 = match(rownames(MetNormal), methyl.db$ID)
normalmet = methyl.db$'Symbol' [m3]
nmethyl_list = data.frame(normalmet)
alsokeep = nmethyl_list[,1] != ""
row.names(MetNormal) = nmethyl_list[1:length(rownames(MetNormal)),]
MetNormal = as.matrix(MetNormal, na.rm = TRUE)
MetNormal = MetNormal[alsokeep,]



which.NA <- function(x) {
  sds = apply(x, 1, sd)
  w = which(is.na(sds))
  names(w)
}

# returns a list containing matrices X1, X2, X3, keeping only the common rows
Make.Common.Rows <- function(X1, X2, X3, rm = NULL) {
  
  common = intersect(rownames(X1), rownames(X2))
  common = intersect(common, rownames(X3))

  if (!is.null(rm)) {
    common = setdiff(common, rm)
  }
  
  m1 = match(common, rownames(X1))
  m2 = match(common, rownames(X2))
  m3 = match(common, rownames(X3))
  
  X1 = X1[m1,]
  X2 = X2[m2,]
  X3 = X3[m3,]
  
  l = list(X1 = X1, X2 = X2, X3 = X3)
  return (l)
}


do.not.keep = c(which.NA(MetNormal), which.NA(MetCancer), which.NA(ExpCancer))

MM = Make.Common.Rows(MetCancer, MetNormal, ExpCancer, rm = do.not.keep)

MetCancerCommon = MM[[1]]
MetNormalCommon = MM[[2]]
ExpCancerCommon = MM[[3]]




MethylMixResults = MethylMix(MetCancerCommon, MetNormalCommon, ExpCancerCommon)
