#when plot.best in findbestprobe is set to FALSE
library(GEOquery)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
load("data/TCGA.RData")
source("findbestprobe.R")

geneChoices = unique(as.character(GPL8490$Symbol))
geneChoices2 = unique(as.character(GPL9183$Symbol))
geneChoices = sort(unique(c(geneChoices, as.character(rownames(TCGA.tumor)), geneChoices2)))
geneChoices = geneChoices[geneChoices!=""]

i = 1
FClistx1 = NULL 
FDRlistx1 = NULL
FClistx2 = NULL 
FDRlistx2 = NULL
FClistx3 = NULL 
FDRlistx3 = NULL
FClistx4 = NULL 
FDRlistx4 = NULL
genesumlist = NULL
genesum = NULL




calc.score <-function(x) {

  if (is.null(x) == TRUE){
    return (0)
  } else if (is.na(x$FC)  | is.na(x$FDR)) {
    return (0)
  }
  
  if (x$FC > 1 & x$FDR <0.1) {
    return (1)
  }
  if (x$FC < 1 & x$FDR <0.1) {
    return (-1)
  }
  return (0)
}


for (i in 1:length(geneChoices)){ 
  genesum = 0
  x1 = findbestprobe(geneChoices[i], GPL8490, GSE37817.methyl, GSE37817.methyl.tumor, "KRIBAB")
  x2 = findbestprobe(geneChoices[i], GPL8490, GSE33510.meth, GSE33510.meth.tumor_normal, "LU")
  x3 = findbestprobe(geneChoices[i], GPL9183, GSE28094.meth, GSE28094.methyl.tumor, "IUOPA")
  x4 = evaluate.paired(geneChoices[i], TCGA.tumor, TCGA.normal, "TCGA")
  
  genesum = calc.score(x1) + calc.score(x2) + calc.score(x3) + calc.score(x4)
  
  if (is.null(x1) == TRUE || is.na(x1) == TRUE){
    x1$FC = "no value"
    x1$FDR = "no value"
  }

  if (is.null(x2) == TRUE || is.na(x2) == TRUE){
    x2$FC = "no value"
    x2$FDR = "no value"
  }
  if (is.null(x3) == TRUE || is.na(x3) == TRUE){
    x3$FC = "no value"
    x3$FDR = "no value"
  }
  if (is.null(x4) == TRUE || is.na(x4) == TRUE){
    x4$FC = "no value"
    x4$FDR = "no value"
  }
  
  
  
  FClistx1 = c(FClistx1, x1$FC)
  FDRlistx1 = c(FDRlistx1, x1$FDR)
  FClistx2 = c(FClistx2, x2$FC)
  FDRlistx2 = c(FDRlistx2, x2$FDR)
  FClistx3 = c(FClistx3, x3$FC)
  FDRlistx3 = c(FDRlistx3, x3$FDR)
  FClistx4 = c(FClistx4, x4$FC)
  FDRlistx4 = c(FDRlistx4, x4$FDR)
  genesumlist = c(genesumlist, genesum)
}

FDRlistx1 = p.adjust(FDRlistx1, method = "fdr")
FDRlistx2 = p.adjust(FDRlistx2, method = "fdr")
FDRlistx3 = p.adjust(FDRlistx3, method = "fdr")
FDRlistx4 = p.adjust(FDRlistx4, method = "fdr")


methyllist = data.frame(Genedataset1.FC = FClistx1, Genedataset1.FDR = FDRlistx1,
                        Genedataset2.FC = FClistx2, Genedataset2.FDR = FDRlistx2,
                        Genedataset3.FC = FClistx3, Genedataset3.FDR = FDRlistx3,
                        Genedataset4.FC = FClistx4, Genedataset4.FDR = FDRlistx4,
                        Count = genesumlist, row.names = geneChoices)


write.csv(methyllist, file = "methyllist2.csv")


