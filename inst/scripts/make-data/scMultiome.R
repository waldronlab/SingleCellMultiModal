library(MultiAssayExperiment)

ddir <- "~/data/scmm/pbmc_10x"

pbmc <- readRDS(
    file.path(ddir, "mae.rds")
)
vdir <- file.path(ddir, paste0("v", "1.0.0"))
setwd(vdir)

## save colData and sampleMap
pbmc_colData <- colData(pbmc)
save(pbmc_colData, file = "pbmc_colData.rda")

pbmc_sampleMap <- sampleMap(pbmc)
save(pbmc_sampleMap, file = "pbmc_sampleMap.rda")

Matrix::writeMM(assay(pbmc[[1]]), "pbmc_rna.mtx")
R.utils::gzip(filename = "pbmc_rna.mtx", destname = "pbmc_rna.mtx.gz")
file.remove("pbmc_rna.mtx")

stopifnot(file.exists("pbmc_rna.mtx.gz"))

rna_mtx <- .read_mtx("pbmc_rna.mtx.gz")

## save H5 file and SCE shell
HDF5Array::saveHDF5SummarizedExperiment(pbmc[[1]], dir = "pbmc_rna",
    prefix = "pbmc_rna_", as.sparse = TRUE)

## load SCE shell
rna_sce <- readRDS("./pbmc_rna/pbmc_rna_se.rds")
## replace assay with MTX assay
pbmc_rna_mtx_obj <- BiocGenerics::replaceSlots(
    rna_sce, assays = Assays(SimpleList(counts = rna_mtx))
)

pbmc_rna_h5_obj <-
    HDF5Array::loadHDF5SummarizedExperiment("pbmc_rna", "pbmc_rna_")

Matrix::writeMM(assay(pbmc[[2]]), "pbmc_atac.mtx")
R.utils::gzip(filename = "pbmc_atac.mtx", destname = "pbmc_atac.mtx.gz")
pbmc_atac_mtx <- "pbmc_atac.mtx.gz"
stopifnot(file.exists(pbmc_atac_mtx))
atac_mtx <- .read_mtx("pbmc_atac.mtx.gz")
## save H5 file and SCE shell
HDF5Array::saveHDF5SummarizedExperiment(pbmc[[2]], dir = "pbmc_atac",
    prefix = "pbmc_atac_", as.sparse = TRUE)

## load SCE shell
atac_sce <- readRDS("./pbmc_atac/pbmc_atac_se.rds")
## replace assay with MTX assay
pbmc_atac_mtx_obj <- BiocGenerics::replaceSlots(
    atac_sce, assays = Assays(SimpleList(counts = atac_mtx))
)

## load H5 object
pbmc_atac_h5_obj <-
    HDF5Array::loadHDF5SummarizedExperiment("pbmc_atac", "pbmc_atac_")
