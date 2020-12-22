
## Make the data to distribute to ExperimentHub

## Required packages
library(HDF5Array)
library(BiocFileCache)

## Steps:
## 1. Retrieve the scRNASeq matrices (n=2) from NCBI
## 2. Read the count matrices
## 3. Store the count matrices in an H5 file
## 4. Retrieve the SCP matrix and annotation from Google Drive
## 5. Store the SCP data in Rda files

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
## 3. Store the count matrices in an H5 file
## ------------------------------------------------------- ##
                                 
h5file <- "../.localdata/SingleCellMultiModal/SCoPE2/v1.0.0/scRNAseq_counts1+2.h5"
writeHDF5Array(m1,
               filepath = h5file,
               name = "rna1",
               with.dimnames = TRUE)
writeHDF5Array(m2,
               filepath = h5file,
               name = "rna2",
               with.dimnames = TRUE)

## ------------------------------------------------------- ##
## 4. Retrieve the SCP matrix and annotation from Google Drive
## ------------------------------------------------------- ##

## Download the protein data provided by the Slavov lab
## https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing
scp_exprs <- read.csv("../.localdata/SCP/specht2019/v3/Proteins-processed.csv",
                      row.names = 1)
scp_exprs <- scp_exprs[, colnames(scp_exprs) != "protein"]
## Download the protein data provided by the Slavov lab
## https://drive.google.com/file/d/16vf6rjIsk-oK9naAH6BQnCFrlWnYtJsS/view?usp=sharing
scp_annot <- read.csv("../.localdata/SCP/specht2019/v3/Cells.csv",
                      row.names = 1)
scp_annot <- t(scp_annot)
scp_annot <- DataFrame(scp_annot)

## ------------------------------------------------------- ##
## 5. Store the SCP data in Rda files
## ------------------------------------------------------- ##

scp_exprs_file <- "../.localdata/SingleCellMultiModal/SCoPE2/v1.0.0/SCoPE2_protein_exprs.Rda"
save(scp_exprs, 
     file = scp_exprs_file)
scp_annot_file <- "../.localdata/SingleCellMultiModal/SCoPE2/v1.0.0/SCoPE2_protein_annot.Rda"
save(scp_annot, 
     file = scp_annot_file)

## ------------------------------------------------------- ##
## Conclusion
## ------------------------------------------------------- ##

## These files should be sent to ExperimentHub
h5file
scp_exprs_file
scp_annot_file
