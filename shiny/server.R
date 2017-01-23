load("C:/Users/Windows User/Desktop/BCBET-M_edited/GSE37817.RData")
y = rep(NA, length(GSE37817.p.expr$title))
y[grep("Control", GSE37817.p.expr$title)] = "Control"
y[grep("Primary", GSE37817.p.expr$title)] = "BladderCancer"
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
