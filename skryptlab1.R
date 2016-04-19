source("https://bioconductor.org/biocLite.R")
biocLite("R.matlab")
library(R.matlab)
biocLite("hgu95av2.db")
library(hgu95av2.db)
genes_list = readMat("genes_list.mat")
genes_list = unlist(genes_list)
x <- hgu95av2ENTREZID
genes_entrez <- as.list(x[genes_list[1:20]])

setGeneric(konwersja,)

