# get data from cloudstor
# https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
## ./output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
if (FALSE) {
    library(MultiAssayExperiment)

    ddir <- "~/data/scmm/mouse_gastrulation"

    if (!dir.exists(ddir))
        dir.create(ddir, recursive = TRUE)

#   old
#   "scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds"
    scnmt <- readRDS(
        file.path(ddir, "allcells",
            "scnmtseq_gastrulation_mae_AllCells.rds"
        )
    )

    exportClass(scnmt, ddir, fmt = "csv")
}

if (FALSE) {
    library(MultiAssayExperiment)

    ddir <- "~/data/scmm/pbmc_10x"

    pbmc <- readRDS(
        file.path(ddir, "mae.rds")
    )
    ## TODO: save sampleMap and colData

    vdir <- file.path(ddir, paste0("v", "1.0.0"))
    setwd(vdir)

    Matrix::writeMM(assay(pbmc[[1]]), "pbmc_rna.mtx")
    R.utils::gzip(filename = "pbmc_rna.mtx", destname = "pbmc_rna.mtx.gz")
    file.remove("pbmc_rna.mtx")

    stopifnot(file.exists("pbmc_rna.mtx.gz"))

    rna_mtx <- HCAMatrixBrowser:::.read_mtx("pbmc_rna.mtx.gz")

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
    atac_mtx <- HCAMatrixBrowser:::.read_mtx("pbmc_atac.mtx.gz")
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
}

# convert .csv files to .rda matrices
.convertData <- function(
    directory = "~/data/scmm/",
    dataDir = "mouse_gastrulation",
    pattern = ".csv")
{
    location <- file.path(directory, dataDir)
    csvs <- list.files(location, pattern = pattern, full.names = TRUE,
        recursive = FALSE)
    invisible(
        lapply(csvs, function(csvfile) {
            objname <- gsub(pattern, "", basename(csvfile))
            readin <- as.data.frame(readr::read_csv(csvfile))
            rnames <- readin[[1L]]

            if (!objname %in% c("scnmt_colData", "scnmt_sampleMap"))
                readin <- data.matrix(readin[, -1])
            else if (identical(objname, "scnmt_colData"))
                names(readin)[1] <- "cellID"
            else
                readin <- readin[, -1]

            if (!objname %in% "scnmt_sampleMap")
                rownames(readin) <- rnames

            assign(objname, readin)
            rdafile <- gsub("csv", "rda", csvfile)
            save(list = objname, file = rdafile)
        })
    )
}
