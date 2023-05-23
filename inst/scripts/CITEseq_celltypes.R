library(devtools)
load_all()
mae <- CITEseq("cord_blood", dry.run=FALSE)

## Detecting MOUSE/HUMAN cells
rna <- assays(mae)[["scRNAseq"]]
hrna <- colSums(rna[grep("^HUMAN", rownames(rna)),])
mrna <- colSums(rna[grep("^MOUSE", rownames(rna)),])
mate <- cbind(hrna, mrna)
plot(log1p(hrna+1), log1p(mrna+1), xlab="hrna", ylab="mrna")
## Using kmeans for detecting 2 bigger clusters of human and mouse cells
set.seed(666)
km <- kmeans(cbind(log(hrna+1), log(mrna+1)), centers=2)
plot(log(hrna+1), log(mrna+1), xlab="hrna", ylab="mrna", col=km$cluster)

## computing distance+hclust on human/mouse cells cluster for detecting mixed cells
mat <- cbind(log(hrna+1), log(mrna+1))

## human cells
d <- dist(mat[km$cluster==2,])
hc <- hclust(d, method="single")
cl <- cutree(hc, k=2)
plot(mat[km$cluster==2,], col=cl)
hbc <- names(cl)[cl==1]

# cd <- colData(mae)
load("cord_blood/v1.0.0/coldata_scRNAseq.rda")
cd <- coldata_scRNAseq

cd$species <- NA
cd$species[which(rownames(cd) %in% hbc)] <- "HUMAN"
cd$species[which(rownames(cd) %in% names(cl)[cl==2])] <- "MIXED"
cd$species[which(rownames(cd) %in% names(km$cluster)[km$cluster==1])] <- "MOUSE"
table(cd$specie)

##### Annotating cell types
adtclrgeo <- as.matrix(read.csv("~/Downloads/GSE100866_CBMC_8K_13AB_10X-ADT_clr-transformed.csv", row.names=1))
## add this assay to the ADTs 
cd$celltype <- NA
cd$markers <- NA

cdct <- cd[!cd$discard,]
cdct <- cdct[cdct$species=="HUMAN",]

out.cd19.cd3 <- getCellGroups(adtclrgeo, adt1="CD19", adt2="CD3", th1=0.9, th2=0.6)


cdct <- addCTLabels(cdct, out.cd19.cd3, "CD19-/CD3+", "T-cells")
cdct <- addCTLabels(cdct, out.cd19.cd3, "CD19+/CD3-", "B-cells")

out.cd11.cd14 <- getCellGroups(adtclrgeo, adt1="CD11c", adt2="CD14", th1=0.4, th2=0.55)
cdct <- addCTLabels(cdct, out.cd11.cd14, "CD11c+/CD14+", "Monocytes CD14+")
table(cdct$celltype)

out.cd11.cd16 <- getCellGroups(adtclrgeo, adt1="CD11c", adt2="CD16", th1=0.4, th2=0.55)
cdct <- addCTLabels(cdct, out.cd11.cd16, "CD11c+/CD16+", "Monocytes CD16+")
table(cdct$celltype)

out.T.cd8.cd4 <- getCellGroups(adtclrgeo[,out.cd19.cd3$`CD19-/CD3+`$bc], adt1="CD8", adt2="CD4", th1=0.9, th2=0.6)
## overwriting because CD4/CD8 T-cells are subgroups of T-cells
cdct <- addCTLabels(cdct, out.T.cd8.cd4, "CD8-/CD4+", "CD4 T-cells", overwrite=TRUE)
cdct <- addCTLabels(cdct, out.T.cd8.cd4, "CD8+/CD4-", "CD8 T-cells", overwrite=TRUE)
# cord_blood_colData_anno <- cdct
# save(cord_blood_colData_anno, file="cord_blood_colData_anno.rda")

## precursors are CD34+ and I took CD56- which seems not expressed from the paper figure
out.cd56.cd34 <- getCellGroups(adtclrgeo, adt1="CD56", adt2="CD34", th1=0.37, th2=0.9)
prebc <- out.cd56.cd34$`CD56-/CD34+`$bc
# idxpre <- which(rownames(cdct) %in% prebc)
# which(prebc %in% rownames(cdct)[!is.na(cdct$celltype)]) ## showing overlap!!!
cdct <- addCTLabels(cdct, out.cd56.cd34, "CD56-/CD34+", "Precursors")

# ## NATURAL KILLERS are CD3-/CD16+ (CD56+ and CD16+)
out.cd16.cd3 <- getCellGroups(adtclrgeo, adt1="CD16", adt2="CD3", th1=0, th2=0.55)
nkcellbc16 <- out.cd16.cd3$`CD16+/CD3-`$bc
idxnk <- which(rownames(cdct) %in% nkcellbc16)
length(idxnk)
sum(nkcellbc16 %in% rownames(cdct)[!is.na(cdct$celltype)]) ## showing overlap!!!

cdctnk <- addCTLabels(cdct, out.cd16.cd3, "CD16+/CD3-", "Natural Killers")

## other markers for NK are CD56+ and CD3-
out.cd56.cd3 <- getCellGroups(adtclrgeo, adt1="CD56", adt2="CD3", th1=0, th2=0)
# nkcellbc56 <- out.cd56.cd3$`CD56+/CD3-`$bc
# idxnk <- which(rownames(cdct) %in% nkcellbc56)
# length(idxnk)
# sum(nkcellbc56 %in% rownames(cdct)[!is.na(cdct$celltype)]) ## showing overlap!!!
cdctnk <- addCTLabels(cdctnk, out.cd56.cd3, "CD56+/CD3-", "Natural Killers")

coldata_scRNAseq <-  cdctnk
save(coldata_scRNAseq, file="cord_blood/v1.0.0/coldata_scRNAseq.rda")
scADT_clrCounts <- adtclrgeo
save(scADT_clrCounts, file="cord_blood/v1.0.0/scADT_clrCounts.rda")


## Building tsne
cnts <- (mae[["scADT"]][, which(colnames(mae[["scADT"]]) %in% rownames(cdctnk))])
adtclrgeoss <- adtclrgeo[,which(colnames(adtclrgeo) %in% rownames(cdctnk))]
adtsce <- SingleCellExperiment(assays=list(counts=cnts, logcounts=adtclrgeoss), 
                               colData=cdctnk)
library(scran)
adtsce <- runPCA(adtsce)
adtsce <- runTSNE(adtsce, dimred="PCA")
plotReducedDim(adtsce, dimred="TSNE", colour_by="celltype")

