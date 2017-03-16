library(GEOquery)
dd = getGEO("GPL8490")
GPL8490 = Table(dd)

save(GPL8490, file = "GPL8490.RData")

