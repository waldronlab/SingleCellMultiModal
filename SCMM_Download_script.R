#	CITE-seq: Large scale simultaneous measuremnt of epitopes and transcriptomes in single cells
library(GEOquery)
## gsm <- getGeo('GSE100866') 
gsm <- getGEO('GSE100866')
#Whole genome expression set
##
library(Biobase)
example(ExpressionSet)

eset <- expresssionSet

write.table(exprs(eset), file = "test.txt")
## 
## read data

myeset <- read.table("test.txt")

library(SummarizedExperiment)
SummarizedExperiment(assays = list(counts = data.matrix(myeset)))



