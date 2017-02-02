setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

shinyServer(function(input, output) {
  output$GenesPlot <- renderPlot({
    column = input$Genes 
    findprobe = GPLmeth$ID[column]
    matching = match(findprobe, rownames(GSE37817.methyl))
    newvector = c(1)
    i = 1
    for (i in 1:length(find)) {
      m = find[i]
      s = split(GSE37817.methyl[matching,], GSE37817.methyl.tumor)
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
    main = paste(column, "expression")
    boxplot(s,
            col = c(input$Color1,input$Color2), main = main, ylab = input$Probes)
  })
  
})
 

