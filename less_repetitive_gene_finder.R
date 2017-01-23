library(GEOquery)
## methylation datasets
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

getset2 = function() {
  if (platform == "GSE33510" || platform == "GSE37817"){
    return (GPLmeth)
  }

else if (platform == "GSE28094"){
  return (GPLmeth2)
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
set2 = getset2()
matching = grep(gene, set2$Symbol) #Match gene to gene name in platform data
findprobe = set2$ID[matching] #Find probe(s) for gene
find = match(findprobe, rownames(set)) #Match probe(s) to row(s) in patient dataset

findbestprobe = function(){
  newvector = c(1)
  i = 1
  for (i in 1:length(find)) {
    m = find[i]
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
