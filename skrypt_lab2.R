load("ExpSet.RData")
library('affy')
library('Biobase')
library('hgu95av2.db')
library('gahgu95av2.db')

setwd("/home/mstolarczyk/Uczelnia/MGR/WSP/lab2/")

##### load datasetA_scans #####
dataScanFile <- file.path(getwd(), "datasetA_scans.txt")
##### create phenoData #####
pData <- read.table( "datasetA_scans.txt", header=TRUE, sep="\t")
summary(pData)

metadata <- data.frame(labelDescription=c("Anotacja","Tkanka","Probka","Nazwa pliku CEL"),row.names=colnames(pData))
phenoData=new("AnnotatedDataFrame",data=pData, varMetadata=metadata)
##### WCZYTANIE DANYCH MIKROMACIERZOWYCH #####
filePaths=paste(pData$scan,'.CEL',sep='') 
data_raw=ReadAffy(filenames=filePaths,phenoData=phenoData) 
data_raw@annotation=paste("ga",data_raw@annotation,sep="")
data_raw@cdfName=paste("GA",data_raw@cdfName,sep="") 


##### NORMALIZACJA RMA #####
RMA=expresso(data_raw,bgcorrect.method="rma",normalize.method="quantiles",
             pmcorrect.method="pmonly",summary.method="medianpolish")
data=exprs(RMA)
colnames(data)=pData$scan

##### FEATURE DATA ####

fData=data.frame(Gene=c(1:nrow(data)),Description=c(1:nrow(data)),row.names=rownames(data)) # przyk?adowe warto?ci

metadata2 <- data.frame(labelDescription=c("Nazwa genu","Opis genu"),row.names=colnames(fData))

featureData=new("AnnotatedDataFrame",data=fData, varMetadata=metadata2) # Annotated Data Frame z datainf

##### OPIS EKSPERYMENTU #####
experimentData=new("MIAME",
                   name="Classification of Human Lung Carcinomas by mRNA Expression Profiling Reveals Distinct Adenocarcinoma Sub-classes",
                   lab="Harvard Medical School",
                   contact="Arindam_Bhattacharjee@dfci.harvard.edu, staunton@genome.wi.mit.edu, Matthew_Meyerson@dfci.harvard.edu, golub@genome.wi.mit.edu",
                   title="Classification of Human Lung Carcinomas by mRNA Expression Profiling Reveals Distinct Adenocarcinoma Sub-classes",
                   abstract="We have generated a molecular taxonomy of lung carcinoma, the leading cause of cancer death in the United States and worldwide. Using oligonucleotide microarrays, we analyzed mRNA expression levels corresponding to 12,600 transcript sequences in 186 lung tumor samples, including 139 adenocarcinomas resected from the lung. Hierarchical and probabilistic clustering of expression data defined distinct sub-classes of lung adenocarcinoma. Among these were tumors with high relative expression of neuroendocrine genes and of type II pneumocyte genes, respectively. Retrospective analysis revealed a less favorable outcome for the adenocarcinomas with neuroendocrine gene expression. The diagnostic potential of expression profiling is emphasized by its ability to discriminate primary lung adenocarcinomas from metastases of extra-pulmonary origin. These results suggest that integration of expression profile data with clinical parameters could aid in diagnosis of lung cancer patients",
                   url="http://www.broadinstitute.org/mpr/lung/")

sampleNames(data_raw@protocolData)=sampleNames(phenoData)

##### STWORZENIE OBIEKTU ExpressionSet #####
rownames(phenoData)=colnames(data)
expset=ExpressionSet(assayData=data,phenoData=phenoData,
                     experimentData=experimentData,
                     annotation=data_raw@annotation,featureData=featureData,protocolData = data_raw@protocolData)
save.image(file = paste("ExpSet",".RData",sep = ""));

