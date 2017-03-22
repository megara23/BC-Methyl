from Bio import Entrez
Entrez.email = "set e-mail here"

f = open("fixThese.csv")

a = f.readlines()
a = [x.strip() for x in a]

gids = [x.split(":")[1] for x in a[1:]]


updateSymbol = []

for id in gids :
  print "retreiving id=" + id + "..."
  handle = Entrez.esummary(db="gene", id=id)
  results = Entrez.read(handle)  # use read instead of parse!
  results = results['DocumentSummarySet']
  results = results['DocumentSummary']
  name = results[0]['Name']
  updateSymbol.append(name)

f = open("updatedGenes.csv", "w")
f.write(a[0] + " updatedSymbol\n")
for i in range(0,len(updateSymbol)) :
  f.write(a[i+1] + " " + updateSymbol[i] + "\n")
f.close()

