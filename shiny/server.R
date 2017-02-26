library(GEOquery)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/TCGA.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

shinyServer(function(input, output) {
  
  output$SummaryPlot <- renderPlot({
    boxplot(c(1:10), main = paste("Summary Plot"))
  })
  
  output$GenesPlot <- renderPlot({
   gene = input$Genes 
   findbestprobe(gene, GPLmeth, GSE37817.methyl, GSE37817.methyl.tumor, "KRIBAB")
  })
  
  output$GenesPlot2 <- renderPlot({
    gene = input$Genes 
    findbestprobe(gene, GPLmeth, GSE33510.meth, GSE33510.meth.tumor_normal, "LU")
  })
  
  output$GenesPlot3 <- renderPlot({
    gene = input$Genes 
    findbestprobe(gene, GPLmeth2, GSE28094.meth, GSE28094.methyl.tumor, "IUOPA")
  })
  
  output$GenesPlot4 <- renderPlot({
    gene = input$Genes 
    evaluate.paired(gene, TCGA.tumor, TCGA.normal, "TCGA Bladder")
  })
  
 })

 