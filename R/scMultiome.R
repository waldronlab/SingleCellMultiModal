scMultiome <-
    function(
        DataType = "pbmc_10x", modes = "*", version = "1.0.0",
        format = c("HDF5", "MTX"), dry.run = TRUE, verbose = TRUE, ...
    )
{
    stopifnot(.isSingleChar(version), .isSingleChar(DataType))
    meta <- list(call = match.call(), version = version)

    if (version != "1.0.0")
        stop("Invalid 'version'; see '?scMultiome' for details.")

    ess_list <- .getResourcesList(prefix = "pbmc_", datatype = DataType,
        modes = modes, version = version, dry.run = dry.run,
        verbose = verbose, ...)

    if (dry.run) { return(ess_list) }

    MultiAssayExperiment(
        experiments = ess_list[["experiments"]],
        colData = ess_list[["colData"]],
        sampleMap = ess_list[["sampleMap"]],
        metadata = meta
    )
}
