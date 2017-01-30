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

findbestprobe = function(gene, GPL, X, Y){

  gene = paste0("^",gene,"$")  
  matching = grep(gene, GPL$Symbol) #Match gene to gene name in platform data
  findprobe = GPL$ID[matching] #Find probe(s) for gene
  find = match(findprobe, rownames(X)) #Match probe(s) to row(s) in patient dataset
  
  newvector = c(1)
  i = 1
  for (i in 1:length(find)) {
    m = find[i]
    s = split(X[m,], Y)
    boxplot(s, main = "Methylation", col = c("purple", "pink"), ylab = "Probe Expression")
    means = lapply(s, mean)
    meanchange = means[[1]] - means[[2]]
    a = s[[1]]
    b = s[[2]]
    z = t.test(a,b)
    p_value= z$p.value
    print(2**meanchange) #Fold Change
    newvector = c(newvector, p_value)
    
    print(newvector)
    }
  newvector = p.adjust(newvector, method = "fdr")
  which.min(newvector) 
  min(newvector)
}


gene = getgene() #Enter gene
findbestprobe(gene, GPLmeth, GSE37817.methyl, GSE37817.methyl.tumor)
findbestprobe(gene, GPLmeth, GSE33510.meth, GSE33510.meth.tumor_normal)
findbestprobe(gene, GPLmeth2, GSE28094.meth, GSE28094.methyl.tumor)

