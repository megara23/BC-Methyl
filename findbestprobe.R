
# finds best probe for a given gene with platform GPL, 
# methylation data X, and phenotypes in Y
# TO DO: needs to handle case where gene is not found
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
    means = lapply(s, mean, na.rm = TRUE)
    meanchange = means[[1]] - means[[2]]
    a = s[[1]]
    b = s[[2]]
    z = t.test(a,b, na.rm = TRUE)
    p_value= z$p.value
    FC = 2**meanchange #Fold Change
    newvector = c(newvector, p_value)
    print(newvector)
  }
  newvector = p.adjust(newvector, method = "fdr")
  which.min(newvector) 
  FDR= newvector[i+1]
  boxplot(s, main = paste(" FC = ", round(FC,2), "FDR = ", round(FDR,4), X), col = c("purple", "pink"), ylab = "Beta Value", names = c("normal", "bladder cancer"))
}


# evaluates differential methylation for paired (e.g., TCGA) data
# function assumes rownames contain the genes
# TO DO: needs to handle case where gene is not found
evaluate.paired <- function(gene, X.tumor, X.normal) {
  m = match(gene, rownames(X.tumor))
  matplot(rbind(X.normal[m,],X.tumor[m,]), 
          ylab = "Methylation (Beta value)", xaxt = "n", 
          type = "b", pch = 19, col = 1, 
          lty = 1, xlim = c(0.9, 2.1))
          axis(1, at = 1:2, labels = c("Normal", "Tumor"))
}