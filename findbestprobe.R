
# finds best probe for a given gene with platform GPL, 
# methylation data X, and phenotypes in Y
# TO DO: needs to handle case where gene is not found
findbestprobe = function(gene, GPL, X, Y, title){

  cat("analyzing, ", title, "\n")
  
  gene = paste0("^",gene,"$")
  matching = grep(gene, GPL$Symbol) #Match gene to gene name in platform data
  i = 1
 if (length(matching) == 0){
   return(NULL)
 }else{
  findprobe = GPL$ID[matching] #Find probe(s) for gene
  find = match(findprobe, rownames(X)) #Match probe(s) to row(s) in patient dataset
  newvector = NULL
  for (i in 1:length(find)) {
    find = match(findprobe, rownames(X)) #Match probe(s) to row(s) in patient dataset
    m = find[i]
    s = split(X[m,], Y)
    s = lapply(s, na.omit)
    a = s[[1]]
    b = s[[2]]
    z = try(t.test(a,b), silent = TRUE)
    if (is(z, "try-error")){
      return(NULL)
    }
    p_value= z$p.value
    newvector = c(newvector, p_value)
    newvector = na.omit(newvector)
    
  }
  newvector = p.adjust(newvector, method = "fdr")
  
  cat("FDR values=")
  print(newvector)
  
  i = which.min(newvector) 
  m_best = find[i]
  s_best = split(X[m_best,], Y)
  s_best = lapply(s_best, na.omit)
  means = lapply(s_best, mean)
  meanchange = means[[1]] - means[[2]]
  FC = 2**meanchange #Fold Change
  FDR= newvector[i]
  cat("Final FDR is ", FDR, "\n")
  
  if (FDR < 0.001 & FDR > 0 ){
    FDR = "< 0.001"
  }
  else if (FDR > 0.001 | FDR == 0.001)
  {
    FDR = round(FDR, 3)
  }
  
  
  file = paste0(title, ".RData")
  
  # the ith probe is the best probe
    #print(s_best)
  save(s_best, i, find, X, m_best, file = file)
    boxplot(s_best, main = paste(title, " FC = ", round(FC,2), "FDR = ", FDR), col = c("purple", "pink"), ylab = "Beta Value", names = c("normal", "bladder cancer"))

  
      }
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