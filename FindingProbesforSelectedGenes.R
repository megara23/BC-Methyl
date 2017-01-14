library(GEOquery)
## methylation datasets
load("C:/Users/Windows User/Desktop/MethylMix/GSE37817.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GSE33510.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GSE28094.RData")
## methylation platforms
load("C:/Users/Windows User/Desktop/MethylMix/GPL8490.RData")
load("C:/Users/Windows User/Desktop/MethylMix/GPL9183.RData")

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
findprobe = GPLmeth$ID[matching] #Find probe(s) for gene
find = match(findprobe, rownames(GSE33510.meth)) #Match probe(s) to row(s) in patient dataset

##Finding the best probe for GSE33510
newvector = c(1)
i = 1
for (i in 1:length(find)) {
  m = find[i]
  s = split(GSE33510.meth[m,], GSE33510.meth.tumor_normal)
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
min(newvector) #Lowest p-value
##Note: remove NaN values 

##Finding the best probe for GSE37817
newvector = c(1)
i=1  
for (i in 1:length(find)) {
  m = find[i]
  s = split(GSE37817.methyl[m,], GSE37817.methyl.tumor)
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
min(newvector) #Lowest p-value

#Since we are working with different platform data, the probes will be different
#Therefore, we must find the probes differently for the GoldenGate Platform
findgene = grep(gene, GPLmeth2$Symbol)
matchgene = GPLmeth2$ID[findgene]
findagain = match(matchgene, rownames(GSE28094.meth))

newvector = c(1)
i = 1
for (i in 1:length(findagain)) {
  m = findagain[i]
  s = split(GSE28094.meth[m,], GSE28094.methyl.tumor)
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
min(newvector) #Lowest p-value
