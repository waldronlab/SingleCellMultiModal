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
# library(scDblFinder)
# scerna <- scDblFinder(rna, k=2, me)
# plot(mat[km$cluster==2,], col=colData(scerna)$scDblFinder.class)
# ## mouse cells
# d <- dist(mat[km$cluster==1,])
# hc <- hclust(d, method="single")
# cl <- cutree(hc, k=9)
# plot(mat[km$cluster==1,], col=cl)
# mbc <- names(cl)[cl==1]
# km1 <- kmeans(mat[km$cluster==1,], centers=3, nstart=20)
# plot(mat[km$cluster==1,], col=km1$cluster)
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
################################ Used Functions

#' addCTLabels
#'
#' @param cd the \code{colData} \code{DataFrame}
#' @param out list data structure returned by \code{getCellGroups}
#' @param outname character indicating the name of the out data structure
#' @param ct character indicating the celltype to assign in the \code{ctcol}
#' @param mkrcol character indicating the cd column to store the markers 
#' indicated by \code{outname} (default is markers)
#' @param ctcol character indicating the column in cd to store the cell type 
#' indicated by \code{ct} (default is celltype)
#' @param overwrite logical indicating if the cell types have to be overwritten 
#' without checking if detected barcodes were already assigned to other celltypes
#' @param verbose logical for having informative messages during the execution
#'
#' @return an updated version of the cd DataFrame 
#' @export
#'
#' @examples
#' TBD
addCTLabels <- function(cd, out, outname, ct, mkrcol="markers", ctcol="celltype",
                        overwrite=FALSE, verbose=TRUE)
{
    ## adds to input cd colData in the mkrcol the markers indicated by outname
    ## and in the ctcol the celltype indicated in ct
    ## the positions for the barcodes (rows in the cd) are taken in position 
    ## outname from the out structure given by function getCellGroups 
    stopifnot(any(c(mkrcol, ctcol) %in% colnames(cd)))
    stopifnot((outname %in% names(out)))
    
    cellbc <- out[[outname]]$bc
    idxc <- which(rownames(cd) %in% cellbc)
    if (length(idxc) !=0)
    {
        if (overwrite)
        {
            if(verbose) message("Blindly overwriting cell types assignments")
            cd[[mkrcol]][idxc] <- outname
            cd[[ctcol]][idxc] <- ct
        } else {
            ## checking if celltypes are already assigned
            idxnona <- which(!is.na(cd[[mkrcol]][idxc]))
            # don't get why ifelse doesn't work
            # idxcnona <- ifelse(length(idxnona)!=0, idxc[-idxnona], idxc)
            if ( length(idxnona)!=0 ) { 
                idxcnona <- idxc[-idxnona]
                if(verbose) message(length(idxnona), " Barcodes already assigned.\n",
                            "Assigning only ", length(idxcnona), " Barcodes...")
            } else { idxcnona <- idxc }
            if (length(idxcnona)!=0) 
            {
                cd[[mkrcol]][idxcnona] <- outname
                cd[[ctcol]][idxcnona] <- ct
            } else {
                if(verbose) message("All selected Barcodes are already assigned\n",
                            "Look at the overwrite argument to handle a more ",
                            "brutal behaviour")
            }
        }
        
    } else {
        warning("No barcodes in cd detected for the selected ", outname, 
                "\nReturning cd as it is...")
    }
    
    return(cd)
}


#' .plotGatingAdt
#' @description
#' 
#' 
#' @param mat 
#' @param adt1 
#' @param adt2 
#' @param th1 
#' @param th2 
#'
#' @return
#' @keywords internal
.plotGatingAdt <- function(mat, adt1="CD19", adt2="CD3", th1=0.2, th2=0)
{
    plot(x=mat[adt1,], y=mat[adt2,], xlab=adt1, ylab=adt2, 
         main=paste0("Gain plot with x-th: ", th1, " y-th: ", th2))
    abline(v=th1, col="red", lty=2)
    abline(h=th2, col="red", lty=2)
    smoothScatter(x=mat[adt1,], y=mat[adt2,], xlab=adt1, ylab=adt2, 
                  main=paste0("Gain plot with x-th: ", th1, " y-th: ", th2))
    
    abline(v=th1, col="red", lty=2)
    abline(h=th2, col="red", lty=2)
}


#' getCellGroups
#' 
#' @description
#' Shows the cells/barcodes in two different plots (scatter and density) 
#' divinding the space in four quadrant indicated by the two thresholds given
#' as input parameters. 
#' The x/y-axis represent respectively the two ADTs given as input.
#' It returns a list of one element for each quadrant, each with barcodes and
#' percentage (see Value section for details).
#'
#' @param mat matrix of counts or clr transformed counts for ADT data in CITEseq
#' @param adt1 character indicating the name of the marker to plot on the x-axis
#' (default is CD19).
#' @param adt2 character indicating the name of the marker to plot on the y-axis
#' (default is CD3).
#' @param th1 numeric indicating the threshold for the marker on the x-axis
#' (default is 0.2).
#' @param th2 numeric indicating the threshold for the marker on the y-axis
#' (default is 0).
#'
#' @return a list of four different element, each one indicating the quarter 
#' where the thresholds divide the plotting space, in eucledian order I, II, 
#' III, IV quadrant, indicating respectively +/+, +/-, -/+, -/- combinations 
#' for the couples of selected ADTs.
#' Each element of the list contains two objects, one with the list of detected 
#' barcodes and one indicating the percentage of barcodes falling into that
#' quadrant.
#' .
#' @details helps to do manual gating for cell type indentification with CITEseq
#' or similar data, providing cell markers. 
#' Once identified two interesting markers for a cell type, the user has to 
#' play with the thresholds to identify the cell populations specified by an 
#' uptake (+) o downtake (-) of the couple of markers (ADTs) previously selected.
#' 
#' @export
#'
#' @examples
#' TBD
getCellGroups <- function(mat, adt1="CD19", adt2="CD3", th1=0.2, th2=0)
{
    stopifnot(any(adt1,adt2) %in% rownames(mat))
    
    plot <- match.arg(plot)
    .plotGatingAdt(mat, adt1, adt2, th1, th2)
    matadt <- mat[c(adt1,adt2),]
    adt1p <- (matadt[adt1,]>th1)
    adt1m <- (matadt[adt1,]<=th1)
    adt2p <- (matadt[adt2,]>th2)
    adt2m <- (matadt[adt2,]<=th2)
    
    
    if (sum(adt1p)+sum(adt1m) != dim(mat)[2]) stop("something went wrong with adt1")
    if (sum(adt2p)+sum(adt2m) != dim(mat)[2]) stop("something went wrong with adt2")
    
    adt12pp <- which(adt1p & adt2p)
    adt12pm <- which(adt1p & adt2m)
    adt12mp <- which(adt1m & adt2p)
    adt12mm <- which(adt1m & adt2m)
    
    l <- list(
        ADT12pp=list(
            bc=colnames(matadt)[adt12pp],
            prc=((length(adt12pp)/dim(matadt)[2])*100)),
        ADT12pm=list(
            bc=colnames(matadt)[adt12pm],
            prc=((length(adt12pm)/dim(matadt)[2])*100)),
        ADT12mp=list(
            bc=colnames(matadt)[adt12mp],
            prc=((length(adt12mp)/dim(matadt)[2])*100)),
        ADT12mm=list(
            bc=colnames(matadt)[adt12mm],
            prc=((length(adt12mm)/dim(matadt)[2])*100))
    )
    names(l) <- c(paste0(adt1,"+/",adt2,"+"),
                  paste0(adt1,"+/",adt2,"-"),
                  paste0(adt1,"-/",adt2,"+"),
                  paste0(adt1,"-/",adt2,"-"))
    
    
    text((min(matadt[adt1,])+0.03), (max(matadt[adt2,])-0.05), paste(round(l[[3]]$prc), "%"))
    text((max(matadt[adt1,])-0.03), (max(matadt[adt2,])-0.05), paste(round(l[[1]]$prc), "%"))
    text((max(matadt[adt1,])-0.03), (min(matadt[adt2,])+0.05), paste(round(l[[2]]$prc), "%"))
    text((min(matadt[adt1,])+0.03), (min(matadt[adt2,])+0.05), paste(round(l[[4]]$prc), "%"))
    return(l)
}


###################################
# 
# library(Seurat)
# # Note that this dataset also contains ~5% of mouse cells, which we can use as negative
# # controls for the protein measurements. For this reason, the gene expression matrix has
# # HUMAN_ or MOUSE_ appended to the beginning of each gene.
# #  To make life a bit easier going forward, we're going to discard all but the top 100 most
# # highly expressed mouse genes, and remove the 'HUMAN_' from the CITE-seq prefix
# seu <- CreateSeuratObject(counts=rna)
# all.equal(colnames(rna), colnames(assays(mae)[["scADT"]]))
# seu[["ADT"]] <-  CreateAssayObject(counts = assays(mae)[["scADT"]])
# DefaultAssay(seu) <- "ADT"
# seu <- NormalizeData(seu, normalization.method = "CLR", margin = 2, assay="ADT")
# # 
# # clr <- seu[["ADT"]]@data
# # plot(clr["CD19",], clr["CD3",])
# # 
# # scaledclr <- scale(seu[["ADT"]]@data)
# # plot(scaledclr["CD19",], scaledclr["CD3",])
# 
# logclradt <- log(seu[["ADT"]]@data)
# min(logclradt)
# plot(logclradt["CD19",], logclradt["CD3",])
# abline(v=0.2)
# abline(h=0)
# 
# out.cd19.cd3 <- getCellGroups(logclradt, adt1="CD19", adt2="CD3",th1=0.2, th2=0)
# 
# tcellbc <- out.cd19.cd3$`CD19-/CD3+`$bc
# length(tcellbc)
# adttcell<-logclradt[,tcellbc]
# dim(adttcell)
# outtcell.cd4.cd8 <- getCellGroups(adttcell, adt1="CD8", adt2="CD4",th1=-0.1, th2=0)
# tcellcd8bc <- outtcell.cd4.cd8$`CD8+/CD4-`$bc
# 
# 
# 
# 
# tcellcd4bc <- outtcell.cd4.cd8$`CD8-/CD4+`$bc
# 
# 
# 
# bcellbc <- tcellbc <- out.cd19.cd3$`CD19è/CD3-`$bc
# 
# out.cd14.cd11 <- getCellGroups(logclradt, adt1="CD11c", adt2="CD14",th1=-0.1, th2=0)
# 
# 
# ###-------
# clrs <- as.matrix(read.csv("~/Downloads/GSE100866_PBMC_vs_flow_10X-ADT_clr-transformed.csv", row.names=1))
# clrs[1:10,1:10]
# out.cd19.cd3 <- getCellGroups(clrs, adt1="CD19", adt2="CD3",th1=1, th2=0.6)
# 
# cnts <- as.matrix(read.csv("~/Downloads/GSE100866_PBMC_vs_flow_10X-ADT_umi.csv", row.names=1))
# cnts[1:10, 1:10]
# 
# ## CBMC
# clrs <- as.matrix(read.csv("/Users/inzirio/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv", row.names=1))
# clrs[1:10, 1:10]
# 
#     

