library(limma)
library(GEOquery)

## methylation datasets
load("C:/Users/Windows User/Desktop/MethylMix/GSE37817.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GSE33510.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GSE28094.RData")

## methylation platforms
load("C:/Users/Windows User/Desktop/MethylMix/GPL8490.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GPL9183.RData")

getprobe = function()
{
  x = readline(prompt = "Enter a probe: ")
  return(as.character(x))
}
probe = getprobe()
  
matching = match(probe, rownames(GSE33510.meth))
s = split(GSE33510.meth[matching,], GSE33510.meth.tumor_normal)
boxplot(s, main = "Methylation", col = c("purple", "pink"), ylab = "Probe Expression")

matching = match(probe, rownames(GSE37817.methyl))
s = split(GSE37817.methyl[matching,], GSE37817.methyl.tumor)
boxplot(s, main = "Methylation", col = c("purple", "pink"), ylab = "Probe Expression")

probe = getprobe() #cg numbers are different for GoldenGate Array

table = Table(GPLmeth2)
findrowname = grep(probe, table$cg_no)
methylprobes = table$'ID'[findrowname]
matching = match(methylprobes, rownames(GSE28094.meth))
s = split(GSE28094.meth[matching,], GSE28094.methyl.tumor)
boxplot(s, main = "Methylation", col = c("purple", "pink"), ylab = "Probe Expression")
