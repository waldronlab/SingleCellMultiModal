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
