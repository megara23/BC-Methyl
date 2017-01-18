library(GEOquery)
## methylation datasets
load("C:/Users/Owner/Desktop/MethylMix/GSE37817.RData")
load("C:/Users/Owner/Desktop/MethylMix/GSE33510.RData")
load("C:/Users/Owner/Desktop/MethylMix/GSE28094.RData")

## methylation platforms
load("C:/Users/Owner/Desktop/MethylMix/GPL8490.RData")
load("C:/Users/Owner/Desktop/MethylMix/GPL9183.RData")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

getgene = function()
{
  x = readline(prompt = "Enter a gene: ")
  return(as.character(x))
}

getplatform = function()
{
  x = readline(prompt = "Enter a dataset: GSE33510, GSE37817, or GSE28094   ")
  return (x)
}

getset = function()
{
  if (platform == "GSE33510"){
    return (GSE33510.meth)
  }
  else if (platform == "GSE37817"){
    return (GSE37817.methyl)
  }
  else if (platform == "GSE28094"){
    return (GSE28094.meth)
  }
  else{
    return ("Not a valid dataset")
  }
}

gettumornorm = function(){
  if (platform == "GSE33510"){
    return (GSE33510.meth.tumor_normal)
  }
  else if (platform == "GSE37817"){
    return (GSE37817.expr.tumor)
  }
  else if (platform == "GSE28094"){
    return (GSE28094.methyl.tumor)
  }
  else{
    return ("Not valid")
  }
}

gene = getgene() #Enter gene
gene = paste0("^",gene,"$") 
platform = getplatform()
set = getset()
matching = grep(gene, GPLmeth$Symbol) #Match gene to gene name in platform data
findprobe = GPLmeth$ID[matching] #Find probe(s) for gene
find = match(findprobe, rownames(GSE33510.meth)) #Match probe(s) to row(s) in patient dataset

#Since we are working with different platform data, the probes will be different
#Therefore, we must find the probes differently for the GoldenGate Platform
findgene = grep(gene, GPLmeth2$Symbol)
matchgene = GPLmeth2$ID[findgene]
findagain = match(matchgene, rownames(GSE28094.meth))

platform_rownames = rownames(set)

pickone = function(){
  if (platform_rownames[1] == "AATK_E63_R"){
    return (findagain)
  }
  else if (platform_rownames[1] == "cg00000292"){
    return(find)
  }
  else {
    return ("error")
  }
}
matched = pickone()

findbestprobe = function(){
  newvector = c(1)
  i = 1
  for (i in 1:length(matched)) {
    m = matched[i]
    s = split(set[m,], gettumornorm())
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
result = findbestprobe()
