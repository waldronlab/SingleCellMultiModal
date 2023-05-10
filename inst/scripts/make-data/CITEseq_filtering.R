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
pb <- CITEseq("peripheral_blood", dry.run=FALSE, DataClass="SingleCellExperiment")
adt <- SingleCellExperiment(assays=list(counts=assays(altExp(pb))[[1]]))
condition <- unlist(lapply(strsplit(colnames(adt), "_"), function(x) x[1]))
colData(adt) <- DataFrame("barcodes"=colnames(adt), "condition"=condition)

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
condition <- unlist(lapply(strsplit(colnames(pb), "_"), function(x) x[1]))
colData(pb) <- DataFrame("barcodes"=colnames(pb), "condition"=condition)
pbcr <- pb[,pb$condition=="CTRL"]
pbcl <- pb[,pb$condition=="CTCL"]
mito <- grep("mt-", tolower(rownames(pb)))
dfcr <- perCellQCMetrics(pbcr, subsets=list(Mito=mito))
dfcl <- perCellQCMetrics(pbcl, subsets=list(Mito=mito))
mito.discardcr <- isOutlier(dfcr$subsets_Mito_percent, type="higher")
summary(mito.discardcr)

mito.discardcl <- isOutlier(dfcl$subsets_Mito_percent, type="higher")
summary(mito.discardcl)


discardcr <- qc.statscr$discard | mito.discardcr
discardcl <- qc.statscl$discard | mito.discardcl

tctrl <- rep(FALSE, dim(pbcl)[2])
length(tctrl)
tctcl <- rep(FALSE, dim(pbcr)[2])
length(tctcl)

discard <- c(discardcr, discardcl)

colData(pb) <- cbind.DataFrame(colData(pb), 
    adt.discard.CTCL=c(qc.statscl$discard, tctcl), 
    mito.discard.CTCL=c(mito.discardcl, tctcl),
    adt.discard.CTRL=c(tctrl, qc.statscr$discard), 
    mito.discard.CTRL=c(tctrl, mito.discardcr),
    discard=discard)

sum(c(qc.statscl$discard, tctcl))
sum(c(mito.discardcl, tctcl))

