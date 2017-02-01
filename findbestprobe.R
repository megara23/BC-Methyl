
# finds best probe for a given gene with platform GPL, 
# methylation data X, and phenotypes in Y
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
