
## Make the data to distribute to ExperimentHub

## Required packages
library(HDF5Array)
library(BiocFileCache)
library(SingleCellExperiment)

## Steps:
## 1. Retrieve the scRNASeq matrices (n=2) from NCBI
## 2. Read the count matrices
## 3. Combine the two RNA count matrices in a SingleCellExperiment
## 4. Retrieve the protein matrix and annotation from Google Drive
## 5. Combine the protein data in a SingleCellExperiment object

## Note that step 3 is optional. I needed to migrate to HDF5 due to 
## memory limitations.

## -------------------------------------- ##
## 1. Retrieve the scRNASeq matrices (n=2) from NCBI
## -------------------------------------- ##

## See also https://bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/UsingHDF5Array.html
bfc <- BiocFileCache()
url1 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4226nnn/GSM4226877/suppl/GSM4226877_rna_data_Bio_Replicate_1.csv.gz"
url2 <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4226nnn/GSM4226878/suppl/GSM4226878_rna_data_Bio_Replicate_2.csv.gz"
bfcrpath(bfc, url1)
bfcrpath(bfc, url2)

## ---------------------- ##
## 2. Read count matrices
## ---------------------- ##

## Batch 1
m1 <- data.table::fread(file = bfcquery(bfc, "GSM4226877")$rpath,
                        sep = ",", header = TRUE)
rn <- m1[[1]]
m1 <- as.matrix(m1[, -1])
rownames(m1) <- rn
colnames(m1) <- paste0(colnames(m1), ".1")
## Batch 2
m2 <- data.table::fread(file = bfcquery(bfc, "GSM4226878")$rpath,
                        sep = ",", header = TRUE)
rn <- m2[[1]]
m2 <- as.matrix(m2[, -1])
colnames(m2) <- paste0(colnames(m2), ".2")
rownames(m2) <- rn

## ------------------------------------------------------- ##
## 3. Combine the two RNA count matrices in a SingleCellExperiment
## ------------------------------------------------------- ##

m1 <- DelayedArray(m1)
m2 <- DelayedArray(m2)
m3 <- cbind(m1, m2) ## This process is delayed until writing
batch <- factor(gsub("^.*[.](\\d)$", "\\1", colnames(m3)))
sce <- SingleCellExperiment(list(counts = m3),
                            colData = DataFrame(Batch = batch))
## The object is rather big and is better stored on disk as an HDF5
saveHDF5SummarizedExperiment(sce, 
                             dir = "../.localdata/SingleCellMultiModal/macrophage_differentiation/v1.0.0/",
                             prefix = "macrophage_rna_", 
                             as.sparse = TRUE)
## Restore some RAM
rm(m1, m2, m3); gc()

## ------------------------------------------------------- ##
## 4. Retrieve the protein matrix and annotation from Google Drive
## ------------------------------------------------------- ##

## Download the protein data provided by the Slavov lab
## https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing
protein_assay <- read.csv("../.localdata/SCP/specht2019/v3/Proteins-processed.csv",
                          row.names = 1)
protein_assay <- protein_assay[, colnames(protein_assay) != "protein"]
protein_assay <- as.matrix(protein_assay)

## Download the protein data provided by the Slavov lab
## https://drive.google.com/file/d/16vf6rjIsk-oK9naAH6BQnCFrlWnYtJsS/view?usp=sharing
protein_colData <- read.csv("../.localdata/SCP/specht2019/v3/Cells.csv",
                            row.names = 1)
protein_colData <- t(protein_colData)
protein_colData <- DataFrame(protein_colData)
## Replace the cell type annotation by more explicit values
protein_colData$celltype <- 
    ifelse(protein_colData$celltype == "sc_m0", "Macrophage", "Monocyte")
## Rename the `raw.file` value by `Batch`
colnames(protein_colData)[5] <- "batch_MS"

## ------------------------------------------------------- ##
## 6. Combine the data in a SingleCellExperiment object
## ------------------------------------------------------- ##

macrophage_protein <- SingleCellExperiment(assay = list(logexprs = protein_assay),
                                           colData = protein_colData)
format(object.size(macrophage_protein), "MB")
## Note the protein data can easily fit in memory. We save it as an Rda
save(macrophage_protein, file = "../.localdata/SingleCellMultiModal/macrophage_differentiation/v1.0.0/macrophage_protein.Rda")

## ------------------------------------------------------- ##
## Conclusion
## ------------------------------------------------------- ##

## These files should be sent to ExperimentHub:
##  mRNA
##  - macrophage_rna_assays.h5
##  - macrophage_rna_se.rds
##  Protein
##  - macrophage_protein.Rda
