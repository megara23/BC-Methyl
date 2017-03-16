# The Methylation_Preprocess data was downloaded from http://firebrowse.org/ on 2/2/2017

# read in data, samples are labeled using standard TCGA bar codes with
# 01 denoting solid tumor
# 06 denoting metastatic tumor
# 11 denoting solid tissue normal
dd = read.delim("BLCA.meth.by_min_expr_corr.data.txt")
rownames(dd) = dd$Hybridization.REF
dd = dd[-1,-(1:4)]

# get tumor samples #
g = grep("01$", colnames(dd))
X.tumor = dd[,g]

# get normal samples
g = grep("11$", colnames(dd))
X.normal = dd[,g]

# remove tumor/normal labels 
colnames(X.normal) = gsub(".11$", "", colnames(X.normal))
colnames(X.tumor) = gsub(".01$", "", colnames(X.tumor))

# keep only paired samples
keep = intersect(colnames(X.normal), colnames(X.tumor))
m1 = match(keep, colnames(X.normal))
m2 = match(keep, colnames(X.tumor))
X.normal = X.normal[,m1]
X.tumor = X.tumor[,m2]

# save the results
TCGA.normal = apply(X.normal, 1:2, as.double)
TCGA.tumor = apply(X.tumor, 1:2, as.double)

save(TCGA.normal, TCGA.tumor, file = "TCGA.RData")

