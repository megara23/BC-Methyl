library(shiny)
library(GEOquery)
GSE37817 = getGEO("GSE37817")
GSE37817.expr = exprs(GSE37817[[1]])
choices = c("red", "green", "blue", "purple", "orange")
shinyUI(
  fluidPage(
    titlePanel("Probe Expression"),
    sidebarLayout(      
      sidebarPanel(
        selectInput("Probes", "Probe:", 
                    choices=rownames(GSE37817.expr)[1:50]),
        hr(),
        
        selectInput("Color1", "1st Color:", 
                    choices=choices),
        selectInput("Color2", "2nd Color:", 
                    choices=choices)
         ),
      mainPanel(
        plotOutput("ProbesPlot")  
      )
  )
)
)
