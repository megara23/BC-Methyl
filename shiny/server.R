library(GEOquery)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

shinyServer(function(input, output) {
  
  output$GenesPlot <- renderPlot({
   gene = input$Genes 
   findbestprobe(gene, GPLmeth, GSE37817.methyl, GSE37817.methyl.tumor)
  })
  
  output$GenesPlot2 <- renderPlot({
    gene = input$Genes2 
    findbestprobe(gene, GPLmeth, GSE37817.methyl, GSE37817.methyl.tumor)
  })
  
})
 