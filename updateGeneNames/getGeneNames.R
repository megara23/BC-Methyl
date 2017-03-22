## GPL8490 gene symbols are incorrect for 'date-like' genes
## this script gets the incorrect gene names

load("data/GPL8490.RData")
keep = as.character(GPL8490$Symbol)<"A" & GPL8490$Symbol!=""
symbols = as.character(GPL8490$Symbol[keep])
gids = as.character(GPL8490$Gene_ID)[keep]
ids = as.character(GPL8490$ID)[keep]

fixThese <- data.frame(id = ids, currentSymbol = symbols, gid = gids)
write.table(fixThese, file = "fixThese.csv", row.names = FALSE, quote = FALSE)
