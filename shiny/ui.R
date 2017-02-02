choices = c("red", "green", "blue", "purple", "orange")
shinyUI(
  fluidPage(
    titlePanel("Methylation Expression"),
    sidebarLayout(      
      sidebarPanel(
        selectInput("Genes", "Gene:", 
                    choices=GPLmeth$Symbol[1:5]),
        hr(),
        
        selectInput("Color1", "1st Color:", 
                    choices=choices),
        selectInput("Color2", "2nd Color:", 
                    choices=choices)
      ),
      mainPanel(
        plotOutput("GenesPlot")  
      )
    )
  )
)
