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

