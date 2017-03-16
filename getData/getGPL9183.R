library(GEOquery)
dd = getGEO("GPL9183")
GPL9183 = Table(dd)

save(GPL9183, file = "GPL9183.RData")

