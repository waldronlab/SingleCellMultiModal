.loadHDF5 <- function(ehub, filepaths, verbose) {
    matchres <- grepl("_assays\\.[Hh]5|_se\\.[Rr][Dd][Ss]", filepaths)
    filepaths <- filepaths[matchres]
    fact <- .removeExt(filepaths)
    fact <- gsub("_se|_assays", "", fact)
    h5list <- split(filepaths, fact)
    names(h5list) <- fact
    lapply(h5list, function(h5file, fn) {
        if (verbose)
            message("Working on: ", paste(fn, collapse = ",\n "))
        se_h5 <- grep("_se", h5file, value = TRUE)
        se_obj <- query(ehub, se_h5)[[1L]]
        h5data <- grep("_assay", h5file, value = TRUE, ignore.case = TRUE)
        h5fileloc <- query(ehub, h5data)[[1L]]
        h5array <- HDF5Array::HDF5Array(h5fileloc, "assay001", as.sparse = TRUE)
        SummarizedExperiment::`assays<-`(
            x = se_obj, withDimnames = FALSE,
            value = list(counts = h5array)
        )
    }, fn = names(h5list))
}

.loadMTX <- function(ehub, filepaths, verbose) {
    matchres <-
        grepl("\\.[Mm][Tt][Xx]\\.[Gg][Zz]$|_se\\.[Rr][Dd][Ss]$", filepaths)
    filepaths <- filepaths[matchres]
    fact <- .removeExt(filepaths)
    fact <- gsub("_se", "", fact)
    mtxlist <- split(filepaths, fact)
    names(mtxlist) <- fact
    lapply(mtxlist, function(mtxfile, fn) {
        if (verbose)
            message("Working on: ", paste(fn, collapse = ",\n "))
        se_mtx <- grep("_se", mtxfile, value = TRUE)
        mtxdata <- grep("mtx", mtxfile, value = TRUE, ignore.case = TRUE)
        se <- query(ehub, se_mtx)[[1L]]
        mtxfile <- query(ehub, mtxdata)[[1L]]
        mtxf <- HCAMatrixBrowser:::.read_mtx(mtxfile)

        BiocGenerics::replaceSlots(
            se, assays = Assays(SimpleList(counts = mtxf))
        )
    }, fn = names(mtxlist))
}

#' Single-cell Multiome ATAC + Gene Expression
#'
#' @description scMultiome currently allows users to download 10K Peripheral
#' Blood Mononuclear Cells provided by
#' [10x Genomics website](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets)
#' (`DataType = "pbmc_10x").
#' This technology enables simultaneous profiling of the transcriptome (using
#' 3’ gene expression) and epigenome (using ATAC-seq) from single cells to
#' deepen your understanding of how genes are expressed and regulated across
#' different cell types. Data prepared by Ricard Argelaguet.
#'
#' @details Users are able to choose from either an `HDF5` or `MTX` file format
#'     as the internal data representation. Choosing `HDF5` gives users a
#'     sparse `HDF5Array` class object. The `MTX` (Matrix Market) format allows
#'     users to load a sparse `dgCMatrix` representation.
#'
#' @inheritParams scNMT
#'
#' @param format Either HDF5 or MTX data formats (default HDF5)
#'
#' @md
#'
#' @examples
#'
#' scMultiOme(DataType = "pbmc_10x", modes = "*", dry.run = TRUE)
#'
#' @export
scMultiome <-
    function(
        DataType = "pbmc_10x", modes = "*", version = "1.0.0",
        format = c("HDF5", "MTX"), dry.run = TRUE, verbose = TRUE, ...
    )
{
    stopifnot(.isSingleChar(version), .isSingleChar(DataType))

    format <- match.arg(format)
    meta <- list(call = match.call(), version = version)

    if (version != "1.0.0")
        stop("Invalid 'version'; see '?scMultiome' for details.")

    ess_list <- .getResourcesList(prefix = "pbmc_", datatype = DataType,
        modes = modes, version = version, dry.run = dry.run,
        verbose = verbose, format = format, ...)

    if (dry.run) { return(ess_list) }

    MultiAssayExperiment(
        experiments = ess_list[["experiments"]],
        colData = ess_list[["colData"]],
        sampleMap = ess_list[["sampleMap"]],
        metadata = meta
    )
}
