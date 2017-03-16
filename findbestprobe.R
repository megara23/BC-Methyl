# function to format the FDR (or p-value)
formatFDR <- function(FDR) {
  if (FDR < 0.001 & FDR > 0 ){
    FDR = "<0.001"
  }
  else if (FDR > 0.001 | FDR == 0.001)
  {
    FDR = round(FDR, 3)
    FDR = paste0("=", FDR)
  }
  return (FDR)
}
# reduces margins
setMargins <-function() {
  # sets bottom, left, top, and right margins
  par(mar = c(3, 3.7, 3, .5) + 0.1)
}
# function to create plot when gene is not found
plotGeneNotFound <-function(title) {
  setMargins()
  plot(c(0,1), type = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "",
       main = paste0("\n",title))
  #legend("center", "Gene not found", bty = "n", cex = 1.5)  
  text(1.5, .5, "Methylation data\nunavailable", cex = 1.5)
}


# finds best probe for a given gene with platform GPL, 
# methylation data X, and phenotypes in Y
# TO DO: needs to handle case where gene is not found
findbestprobe = function(gene, GPL, X, Y, title){

  gene = paste0("^",gene,"$")
  matching = grep(gene, GPL$Symbol) #Match gene to gene name in platform data
  i = 1
 if (length(matching) == 0){
   plotGeneNotFound(title)
   return(NULL)
 }else{
  findprobe = GPL$ID[matching] #Find probe(s) for gene
  find = match(findprobe, rownames(X)) #Match probe(s) to row(s) in patient dataset
  length = length(find)
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
      return( plotGeneNotFound(title))
    }
    p_value= z$p.value
    newvector = c(newvector, p_value)
    #newvector = na.omit(newvector)
  }
  
  numFound = sum(!is.na(newvector))
  if (numFound == 0) {
    plotGeneNotFound(title)
    return()
  }
  newvector = p.adjust(newvector, method = "fdr")

  i = which.min(newvector) 
  m_best = find[i]
  s_best = split(X[m_best,], Y)
  s_best = lapply(s_best, na.omit)

  # calculate FC, if NaN, then we have 0/0
  means = lapply(s_best, mean)
  meanchange = means[[2]] / means[[1]]
  FC = meanchange #Fold Change
  if (is.nan(FC)) FC = 1
  
  FDR= newvector[i]
  strFDR = formatFDR(FDR) 
  if (numFound==1){
    pvaltitle = "(P"}
  else{
    pvaltitle = "(FDR"
  }
  # the ith probe is the best probe
  if (!is.null(s_best)) {
    setMargins()
    boxplot(s_best, main = paste(title, "\nFC = ", round(FC,2), pvaltitle, strFDR, ")"), col = c("purple", "pink"), ylab = "Beta Value", names = c("Normal", "Tumor"))
  }
      }
  cat("Final FDR is ", FDR, "\n")
 
  # the ith probe is the best probe
  if (!is.null(s_best)) {
    setMargins()
  cat("Final FDR is ", FDR, "\n")
 
  strFDR = formatFDR(FDR) 
  
  # the ith probe is the best probe
  if (!is.null(s_best)) {
    setMargins()
    boxplot(s_best, main = paste0(title, "\nFC=", round(FC,2), " (FDR", strFDR, ")"), col = c("purple", "pink"), ylab = "Beta Value", names = c("Normal", "Tumor"),
            ylim = c(0,1))
    results = list(FC = FC, FDR = FDR)
    return (results)
  }
 }
}

# evaluates differential methylation for paired (e.g., TCGA) data
# function assumes rownames contain the genes
# TO DO: needs to handle case where gene is not found
evaluate.paired <- function(gene, X.tumor, X.normal, title) {
  m = match(gene, rownames(X.tumor))
  if (is.na(m)){
   plotGeneNotFound(title)
   return(NULL)
  }else{

  t = t.test(X.normal[m,], X.tumor[m,], paired = TRUE) 

  strP = formatFDR(t$p.value)
  FC = mean(X.tumor[m,] / X.normal[m,])
  if (is.nan(FC)) FC = 1
  
  title = paste0(title, "\nFC=", round(FC,2), " (P", strP, ")")
  setMargins()
  matplot(rbind(X.normal[m,],X.tumor[m,]), 
          main = title, ylab = "Beta value", xaxt = "n", 
          type = "b", pch = 19, col = 1, 
          lty = 1, xlim = c(0.9, 2.1), ylim = c(0, 1))
          axis(1, at = 1:2, labels = c("Normal", "Tumor"))
  }
  results = list(FC = FC, FDR = t$p.value)
  return (results)
}
