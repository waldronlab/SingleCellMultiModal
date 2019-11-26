#	CITE-seq: Large scale simultaneous measuremnt of epitopes and transcriptomes in single cells
library(GEOquery)
# gsm <- getGeo('GSM2695379') 
gsm <- getGEO('GSM2695379')
#ADT: CBMC_8K_13AB_10x info for above script
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

