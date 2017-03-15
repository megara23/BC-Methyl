library(GEOquery)
library(ggplot2)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/TCGA.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")
source("plotSummary.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

shinyServer(function(input, output) {
  
  DM.results <-reactiveValues(RESULTS = list(KRIBAB = NULL, LU = NULL, IUOPA = NULL, TCGA = NULL))
  
  output$SummaryPlot <- renderPlot({
    setMargins()
    plotSummary(DM.results$RESULTS)
    
  })
  
  output$GenesPlot <- renderPlot({
   gene = input$Genes 
   DM.results$RESULTS$KRIBAB = findbestprobe(gene, GPLmeth, GSE37817.methyl, GSE37817.methyl.tumor, "KRIBAB")
  })
  
  output$GenesPlot2 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$LU = findbestprobe(gene, GPLmeth, GSE33510.meth, GSE33510.meth.tumor_normal, "LU")
  })
  
  output$GenesPlot3 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$IUOPA = findbestprobe(gene, GPLmeth2, GSE28094.meth, GSE28094.methyl.tumor, "IUOPA")
  })
  
  output$GenesPlot4 <- renderPlot({
    gene = input$Genes 
    DM.results$RESULTS$TCGA = evaluate.paired(gene, TCGA.tumor, TCGA.normal, "TCGA")
  })
  
 })

 
