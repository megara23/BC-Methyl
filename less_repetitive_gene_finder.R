library(GEOquery)
## methylation datasets

setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")

## methylation platforms
load("data/GPL8490.RData")
load("data/GPL9183.RData")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

getgene = function()
{
  x = readline(prompt = "Enter a gene: ")
  return(as.character(x))
}

gene = getgene() #Enter gene
gene = paste0("^",gene,"$") 
matching = grep(gene, GPLmeth$Symbol) #Match gene to gene name in platform data
matching2 = grep(gene, GPLmeth2$Symbol)
findprobe = GPLmeth$ID[matching] #Find probe(s) for gene
findprobe2 = GPLmeth2$ID[matching]
find = match(findprobe, rownames(GSE33510.meth)) #Match probe(s) to row(s) in patient dataset
find2 = match(findprobe2, rownames(GSE33510.meth))
findbestprobe = function(methyl, set, gettumornorm){
  newvector = c(1)
  i = 1
  for (i in 1:length(methyl)) {
    m = methyl[i]
    s = split(set[m,], gettumornorm)
    boxplot(s, main = "Methylation", col = c("purple", "pink"), ylab = "Probe Expression")
    means = lapply(s, mean)
    meanchange = means[[1]] - means[[2]]
    a = s[[1]]
    b = s[[2]]
    z = t.test(a,b)
    p_value= z$p.value
    print(2**meanchange) #Fold Change
    newvector = append(newvector, p_value, after=length(newvector))
    print(newvector)
    }
  newvector = p.adjust(newvector, method = "fdr")
  which.min(newvector) 
  min(newvector)
}
findbestprobe(find, GSE37817.methyl, GSE37817.methyl.tumor)
findbestprobe(find, GSE33510.meth, GSE33510.meth.tumor_normal)
findbestprobe(find2, GSE28094.meth, GSE28094.methyl.tumor)