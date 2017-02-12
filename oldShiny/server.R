load("data/GSE37817.RData")

shinyServer(function(input, output) {
  output$ProbesPlot <- renderPlot({
    column = input$Probes 
    matching = match(input$Probes, rownames(GSE37817.methyl))
    s = split(GSE37817.methyl[matching,], GSE37817.methyl.tumor)
    main = paste(column, "expression")
    boxplot(s,
            col = c(input$Color1,input$Color2), main = main, ylab = input$Probes)
  })
  
})
