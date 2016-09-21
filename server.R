library(shiny)
library(GEOquery)
GSE37817 = getGEO("GSE37817")
GSE37817.expr = exprs(GSE37817[[1]])

GSE37817.p = pData(GSE37817[[1]])
y = rep(NA, length(GSE37817.p$title))
y[grep("Control", GSE37817.p$title)] = "Control"
y[grep("Primary", GSE37817.p$title)] = "BladderCancer"
disorder = as.character(y)
table(disorder)

shinyServer(function(input, output) {
  output$ProbesPlot <- renderPlot({
    column = input$Probes 
    matching = match(input$Probes, rownames(GSE37817.expr))
    s = split(GSE37817.expr[matching,], disorder)
    main = paste(column, "expression")
    boxplot(s,
            col = c(input$Color1,input$Color2), main = main, ylab = input$Probes)
  })
  
})
