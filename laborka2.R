source("http://bioconductor.org/biocLite.R")
biocLite(c("Biobase", "affy"))
library("affy")

biocLite(c("hgu95av2.db"))
biocLite("gahgu95av2.db")
library("hgu95av2.db")
library("gahgu95av2.db")

#########   ZAD 1
#data<- ReadAffy()
dane=read.AnnotatedDataFrame("trzy.txt",row.names = NULL)

#RAW_DATA=read.affybatch(paste(dane@data$scan, sep="", ".CEL")
RAW_DATA=read.affybatch(dane@data$scan)
fpenodaty=pData(dane)
fenodaty=phenoData(data)

#########   ZAD 2 Normalizacja danych (MAS, RMA, GCRMA)
#########   ZAD 3 Annotacje danych mikromacierzowych (Affymetrix, Ferrari)

RAW_DATA@cdfName=paste("ga", RAW_DATA@cdfName,sep="")
RAW_DATA@annotation=paste("ga", RAW_DATA@annotation, sep="")
znormalizowane=rma(RAW_DATA)

#########   ZAD 4 Tworzenie obiektu typu ExpressionSet 
macierz_ekspresji=exprs(znormalizowane)


#tworzenie feature data
annotation <- "hgu95av2"
feature_data<- RAW_DATA@featureData

new_espr_set <- new("ExpressionSet", exprs = macierz_ekspresji, featureData=feature_data)

