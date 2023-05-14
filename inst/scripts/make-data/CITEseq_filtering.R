library(SingleCellExperiment)
library(DropletUtils)
cb <- CITEseq("cord_blood", dry.run=FALSE, DataClass="SingleCellExperiment")
adt <- SingleCellExperiment(assays=list(counts=assays(altExp(cb))[[1]]))
top.marker <- rownames(adt)[max.col(t(counts(adt)))]
total.count <- colSums(counts(adt))
boxplot(split(log10(total.count), top.marker), ylab="Log-total ADT count", las=2)
adt.counts <- counts(adt)
adt.detected <- colSums(adt.counts > 0)
hist(adt.detected, col='grey', main="", xlab="Number of detected ADTs")

qc.stats <- cleanTagCounts(adt)#, exclusive=c("CD3", "CD19"))
summary(qc.stats$high.ambient) # libraries removed with high ambient contamination

library(scater)
mito <- grep("mt-", tolower(rownames(adt)))
df <- perCellQCMetrics(adt, subsets=list(Mito=mito))
mito.discard <- isOutlier(df$subsets_Mito_percent, type="higher")
summary(mito.discard)

discard <- qc.stats$discard | mito.discard

colData(cb) <- cbind.DataFrame(colData(cb), adt.discard=qc.stats$discard, mito.discard=mito.discard, discard=discard)

scRNAseq_coldata <- colData(cb)
dir.create("cord_blood/v1.0.0/", recursive=TRUE)
save(scRNAseq_coldata, file="cord_blood/v1.0.0/scRNAseq_coldata.rda")

## Alternatively it is possible to indicate two or more ADTs that should 
## be expressed alternatively in a cell.
## for CD3/CD4/CD8 i'm referring to https://tinyurl.com/ys9aawce
## otherwise the OSCA vignette 12.3.3 indicates to use CD3/CD19, but
## this article found that CD3and CD19 could be expressed in a novel cell type
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8694500/
# qc.stats <- cleanTagCounts(adtsce1, exclusive=c("CD4", "CD8"))
# summary(qc.stats$discard) # libraries removed with high ambient contamination

library(SingleCellExperiment)
library(DropletUtils)
mae <- CITEseq("peripheral_blood", dry.run=FALSE)#, DataClass="SingleCellExperiment")
adt <- SingleCellExperiment(assays=list(counts=mae[["scADT"]]))
pb <-  SingleCellExperiment(assays=list(counts=mae[["scRNA"]]))
cn <- colnames(adt)
condition <- unlist(lapply(strsplit(colnames(adt), "_"), function(x) x[1]))
bc <- unlist(lapply(strsplit(colnames(adt), "_"), function(x) x[2]))
colData(adt) <- DataFrame("barcodes"=bc, "condition"=condition)
colnames(adt) <- cn

adt.rm <- adt[-c(3,52),]

adtcr <- adt.rm[, adt.rm$condition=="CTRL"]
adtcl <- adt.rm[, adt.rm$condition=="CTCL"]
top.markercr <- rownames(adtcr)[max.col(t(counts(adtcr)))]
top.markercl <- rownames(adtcl)[max.col(t(counts(adtcl)))]
total.countcr <- colSums(counts(adtcr))
total.countcl <- colSums(counts(adtcl))

boxplot(split(log10(total.countcr), top.markercr), ylab="Log-total ADT CTRL count", las=2) #CD5
boxplot(split(log10(total.countcl), top.markercl), ylab="Log-total ADT CTCL count", las=2) #CD279

adt.countscr <- counts(adtcr)
adt.detectedcr <- colSums(adt.countscr > 0)
hist(adt.detectedcr, col='grey', main="", xlab="Number of detected ADTs CTRL")

adt.countscl <- counts(adtcl)
adt.detectedcl <- colSums(adt.countscl > 0)
hist(adt.detectedcl, col='grey', main="", xlab="Number of detected ADTs CTCL")

qc.statscr <- cleanTagCounts(adtcr)#, exclusive=c("CD3", "CD19"))
summary(qc.statscr$high.ambient) # libraries removed with high ambient contamination

qc.statscl <- cleanTagCounts(adtcl)#, exclusive=c("CD3", "CD19"))
summary(qc.statscl$high.ambient) # libraries removed with high am

library(scater)
cn <- colnames(pb)
condition <- unlist(lapply(strsplit(colnames(pb), "_"), function(x) x[1]))
bc <- unlist(lapply(strsplit(colnames(pb), "_"), function(x) x[2]))
colData(pb) <- DataFrame("barcodes"=bc, "condition"=condition)
colnames(pb) <- cn
pbcr <- pb[,pb$condition=="CTRL"]
pbcl <- pb[,pb$condition=="CTCL"]
mito <- grep("mt-", tolower(rownames(pb)))
dfcr <- perCellQCMetrics(pbcr, subsets=list(Mito=mito))
dfcl <- perCellQCMetrics(pbcl, subsets=list(Mito=mito))
mito.discardcr <- isOutlier(dfcr$subsets_Mito_percent, type="higher")
names(mito.discardcr) <- rownames(dfcr)
summary(mito.discardcr)

mito.discardcl <- isOutlier(dfcl$subsets_Mito_percent, type="higher")
names(mito.discardcl) <- rownames(dfcl)
summary(mito.discardcl)

cd <- colData(mae)
cd
cd$adt.discard_CTRL <- FALSE
cd$adt.discard_CTRL[which(rownames(cd) %in% rownames(qc.statscr)[qc.statscr$discard])] <- TRUE
cd$adt.discard_CTCL <- FALSE
cd$adt.discard_CTCL[which(rownames(cd) %in% rownames(qc.statscl)[qc.statscl$discard])] <- TRUE

cd$mito.discard_CTRL <- FALSE
cd$mito.discard_CTRL[which(rownames(cd) %in% names(mito.discardcr)[mito.discardcr])] <- TRUE
cd$mito.discard_CTCL <- FALSE
cd$mito.discard_CTCL[which(rownames(cd) %in% names(mito.discardcl)[mito.discardcl])] <- TRUE
cd$discard_CTRL <- cd$adt.discard_CTRL | cd$mito.discard_CTRL
cd$discard_CTCL <- cd$adt.discard_CTCL | cd$mito.discard_CTCL
cd$discard <- cd$discard_CTRL | cd$discard_CTCL


scRNAseq_coldata <- cd
dir.create("peripheral_blood/v1.0.0/", recursive=TRUE)
save(scRNAseq_coldata, file="peripheral_blood/v1.0.0/scRNAseq_coldata.rda")






