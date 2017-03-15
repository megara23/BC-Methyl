library(GEOquery)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
load("data/TCGA.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

choices = c("red", "green", "blue", "purple", "orange")
geneChoices = unique(as.character(GPLmeth$Symbol))
geneChoices2 = unique(as.character(GPLmeth2$Symbol))
geneChoices = sort(unique(c(geneChoices, as.character(rownames(TCGA.tumor)), geneChoices2)))
geneChoices = geneChoices[geneChoices!=""]

shinyUI(
  fluidPage(
    titlePanel("Methylation Expression"),
    sidebarLayout(      
      sidebarPanel(width = 2,
        selectInput("Genes", "Select a gene:", 
                    choices=geneChoices)
      ),
      mainPanel(
      fluidRow(
        column(4, plotOutput("SummaryPlot", height = 300)),
        column(4, plotOutput("GenesPlot", height = 300)),
        column(4, plotOutput("GenesPlot2", height = 300))
        ),
      fluidRow(
        column(2),
        column(4, plotOutput("GenesPlot3", height = 300)),
        column(4, plotOutput("GenesPlot4", height = 300)),
        column(2)
        
        
        )
      )
    )
  )
)
