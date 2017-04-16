setwd("C:/Users/Owner/Desktop/BC-Methyl/shiny")
load("data/SC.RData")
load("data/Korea.RData")
load("data/GSE3167.RData")
library(GEOquery)


geneChoices = sort(unique(c(as.character(rownames(res.gse3167)), as.character(rownames(res.korea)), as.character(rownames(res.sc)))))
geneChoices = geneChoices[geneChoices!=""]

maketable = function(x){
i = 1 
d = 0
origlength = length(rownames(x))
for (i in (1:length(geneChoices))){
  match = match(geneChoices[i], rownames(x))
  if (is.na(match) == TRUE){
    d = d+1
    x = rbind(x, c("NA", "NA", "NA"))
    row.names(x)[d+origlength] = geneChoices[i]
   }
}
x = x[-2]
x[order(rownames(x)),]
}

res.gse3167 = maketable(res.gse3167)
res.korea = maketable(res.korea)
res.sc = maketable(res.sc)

calc.score <-function(x, i) {
  if (is.na(x[i,1]) == TRUE | is.na(x[i,2]) ==TRUE ){
    return (0)
  }
  
  if (x[i,1] > 1 & x[i,2] <0.1) {
    return (1)
  }
  if (x[i,1] < -1 & x[i,2] <0.1) {
    return (-1)
  }
  return (0)
}


t = 1
genesumlist = NULL 

for (t in 1:length(geneChoices)){
sum = calc.score(res.gse3167, t) + calc.score(res.korea, t) + calc.score(res.sc, t)
genesumlist = c(genesumlist, sum)
}

genelist = data.frame(gseFC= res.gse3167[,1],gseFDR= res.gse3167[,2],korFC= res.korea[,1],korFDR= res.korea[,2],
                      scFC= res.sc[,1],scFDR= res.sc[,2], Count = genesumlist, row.names=geneChoices)

##
setwd("C:/Users/Owner/Desktop/BC-Methyl")
library(gplots)
source("Find_top_methylated_probes.R")

a = c("MAOA", "TP53", "HRAS", "BRCA2")
b = c("FGFR3", "BRCA1", "CD24", "CD44", "MAOA", "ABC")

# how many genes are common to the methylation and expression datasets? (GD can get this)

setdiff(b,a)


counts = list(DE = a, DM = b)
venn(counts)


