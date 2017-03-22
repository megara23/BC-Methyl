## update Gene Symbols for GPL8490 platform 

load("data/GPL8490.RData")
GPL8490$Symbol = as.character(GPL8490$Symbol)
gene = read.csv("updatedGenes.csv", sep = " ", as.is = TRUE)

for (i in 1:nrow(gene)) {
  keep = GPL8490$Symbol%in% gene$currentSymbol[i]
  GPL8490$Symbol[keep] = gene$updatedSymbol[i]
}

save(GPL8490, file = "data/GPL8490.RData")
