library(GEOquery)
library(ggplot2)
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/TCGA.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")
source("plotSummary.R")

shinyServer(function(input, output) {
  
  DM.results <-reactiveValues(RESULTS = list(KRIBAB = NULL, LU = NULL, IUOPA = NULL, TCGA = NULL))
  
  output$SummaryPlot <- renderPlot({
    setMargins()
    plotSummary(DM.results$RESULTS)
    
  })
  
  output$GenesPlot <- renderPlot({
   gene = input$Genes 
   DM.results$RESULTS$KRIBAB = findbestprobe(gene, GPL8490, GSE37817.methyl, GSE37817.methyl.tumor, "KRIBAB")
  })
  
  output$GenesPlot2 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$LU = findbestprobe(gene, GPL8490, GSE33510.meth, GSE33510.meth.tumor_normal, "LU")
  })
  
  output$GenesPlot3 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$IUOPA = findbestprobe(gene, GPL9183, GSE28094.meth, GSE28094.methyl.tumor, "IUOPA")
  })
  
  output$GenesPlot4 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$TCGA = evaluate.paired(gene, TCGA.tumor, TCGA.normal, "TCGA")
  })
  
 })

 
