library(GEOquery)
## methylation datasets

setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")

## methylation platforms
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

#Dataset names 
GSE37817.methyl = KRIBAB 
GSE33510.meth = LUSUC
GSE28094.meth = IUOPA

getgene = function()
{
  x = readline(prompt = "Enter a gene: ")
  return(as.character(x))
}


gene = getgene() #Enter gene
findbestprobe(gene, GPLmeth, KRIBAB, GSE37817.methyl.tumor)
findbestprobe(gene, GPLmeth, LUSUC, GSE33510.meth.tumor_normal)
findbestprobe(gene, GPLmeth2, IUOPA, GSE28094.methyl.tumor)
  
