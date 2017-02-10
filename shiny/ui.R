setwd("C:/Users/Owner/Desktop/BC-Methyl/")
load("data/GSE37817.RData")
load("data/GPL8490.RData")
load("data/GPL9183.RData")
source("findbestprobe.R")

GPLmeth = Table(GPLmeth) #Illumina HumanMethylation27 Beadchip Array Platform Data
GPLmeth2 = Table(GPLmeth2) #GoldenGate Methylation Array Platform Data 

choices = c("red", "green", "blue", "purple", "orange")
shinyUI(
  fluidPage(
    titlePanel("Methylation Expression"),
    sidebarLayout(      
      sidebarPanel(
        selectInput("Genes", "Gene 1:", 
                    choices=GPLmeth$Symbol[1:5]),
        selectInput("Genes2", "Gene 2:", 
                    choices=GPLmeth$Symbol[1:5]),
        hr(),
        
        selectInput("Color1", "1st Color:", 
                    choices=choices),
        selectInput("Color2", "2nd Color:", 
                    choices=choices)
      ),
      mainPanel(
        plotOutput("GenesPlot"),
        plotOutput("GenesPlot2")
      )
    )
  )
)
