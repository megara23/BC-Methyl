library(GEOquery)
setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GSE33510.RData")
load("data/GSE28094.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

choices = c("red", "green", "blue", "purple", "orange")
geneChoices = as.character(GPLmeth$Symbol)
geneChoices2 = as.character(GPLmeth2$Symbol)
geneChoices = c(geneChoices, geneChoices2, "Gene")
shinyUI(
  fluidPage(
    titlePanel("Methylation Expression"),
    sidebarLayout(      
      sidebarPanel(
        selectInput("Genes", "Gene 1:", 
                    choices=geneChoices)
      ),
      mainPanel(
      fluidRow(
        column(4, plotOutput("GenesPlot")),
        column(4, plotOutput("GenesPlot2")),
        column(4, plotOutput("GenesPlot3"))
        )
      )
    )
  )
)
