.getAnswer <- function(msg, allowed)
{
    if (interactive()) {
        repeat {
            cat(msg)
            answer <- readLines(n = 1)
            if (answer %in% allowed)
                break
        }
        tolower(answer)
    } else {
        "n"
    }
}

.isSingleChar <- function(x) {
    length(x) == 1L && is.character(x) && !is.na(x)
}

.removeExt <- function(fnames) {
    gsub("\\..*$", "", basename(fnames))
}

.modesAvailable <- function(listfiles, prefix) {
    slots <- c("metadata", "colData", "sampleMap")
    modes <- gsub(prefix, "", listfiles, fixed = TRUE)
    modes <- gsub("_assays|_se|_tenx", "", modes)
    modes <- .removeExt(modes)
    unique(sort(modes[!modes %in% slots]))
}

.searchFromInputs <- function(glob, searchFields) {
    regGlob <- glob2rx(unique(glob))
    res <- unlist(lapply(regGlob, function(x) {
        grep(x, searchFields, ignore.case = TRUE, value = TRUE)
        }))
    if (!length(res))
        stop("No matches found, modify search criteria")
    res
}

.conditionToIndex <- function(startVec, testVec, FUN) {
    logmat <- vapply(startVec, FUN, logical(length(testVec)))
    apply(logmat, 1L, any)
}

.queryResources <- function(ExperimentHub, resTable, verbose) {
    fileNames <- stats::setNames(resTable[["RDataPath"]], resTable[["Title"]])
    lapply(fileNames, function(res) {
        if (verbose)
            message("Working on: ", gsub("\\.rda", "", basename(res)))
        # only take the last one for multiple matches
        utils::tail(query(ExperimentHub, res), 1)
    })
}

.getResources <- function(ExperimentHub, resTable, prefix, verbose) {
    infos <- .queryResources(ExperimentHub, resTable, verbose)
    rpath <- vapply(infos, function(x) `$`(x, "rdatapath"), character(1L))

    h5resources <- grepl("\\.[Hh]5$", rpath)
    mtxresources <- grepl("\\.[Mm][Tt][Xx]\\.[Gg][Zz]$", rpath)
    shells <- grepl("se\\.[Rr][Dd][Ss]$", rpath)
    otherres <- !((h5resources | mtxresources) | shells)

    if (any(h5resources))
        matress <- .loadHDF5(ExperimentHub, rpath, verbose)
    else if (any(mtxresources))
        matress <- .loadMTX(ExperimentHub, rpath, verbose)
    else
        matress <- list()

    if (any(otherres)) {
        rest <- lapply(infos[otherres], `[[`, 1L)
        c(rest, matress)
    } else {
        matress
    }
}

.getResourceInfo <- function(ExperimentHub, resTable, prefix, verbose) {
    infos <- .queryResources(ExperimentHub, resTable, verbose)
    resID <- vapply(infos, names, character(1L))
    restab <- AnnotationHub::getInfoOnIds(ExperimentHub, resID)
    restab <-
        restab[, !names(restab) %in% c("fetch_id", "status", "biocversion")]
    sizes <- as.numeric(restab[["file_size"]])
    class(sizes) <- "object_size"
    titleidx <- which(names(restab) == "title")
    restab <- as.data.frame(append(
        restab,
        list(mode = gsub(prefix, "", restab[["title"]]),
            file_size = format(sizes, units = "Mb")),
        titleidx
    ))
    restab[, -c(length(restab), titleidx)]
}

.test_eh <- function(...) {
    tryCatch({
        ExperimentHub(...)
    }, error = function(e) {
        emsg <- conditionMessage(e)
        if (grepl("Timeout", emsg))
            warning("[experimenthub.bioconductor.org] timeout, localHub=TRUE",
                call.=FALSE)
        ExperimentHub(..., localHub = TRUE)
    })
}

.isSingleCharNA <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x)
}

.getResourcesList <-
    function(prefix, datatype, modes, version, format, dry.run, verbose, ...)
{
    modes_file <- system.file("extdata", "metadata.csv",
        package = "SingleCellMultiModal", mustWork = TRUE)

    DataType <- tolower(datatype)
    stopifnot(
        .isSingleCharNA(DataType), .isSingleCharNA(version)
    )

    modes_metadat <- read.csv(modes_file, stringsAsFactors = FALSE)
    if (missing(format))
        notfmt <- "FakeFormatNoMatch"
    else
        notfmt <- switch(format, HDF5 = "MTX", MTX = "HDF5", format)
    filt <- modes_metadat[["DataType"]] == DataType &
        modes_metadat[["SourceVersion"]] == version &
        modes_metadat[["SourceType"]] != notfmt

    modes_metadat <- modes_metadat[filt, , drop = FALSE]
    eh_assays <- modes_metadat[["ResourceName"]]
    modesAvail <- .modesAvailable(eh_assays, prefix)
    resultModes <- .searchFromInputs(modes, modesAvail)
    fileIdx <- .conditionToIndex(
        resultModes, eh_assays, function(x) grepl(x, eh_assays)
    )
    fileMatches <- modes_metadat[fileIdx, c("Title", "DispatchClass", "SourceVersion")]
    eh <- .test_eh(...)

    if (dry.run) {
        return(.getResourceInfo(
            eh, modes_metadat[fileIdx, c("Title", "RDataPath")], prefix, FALSE
        ))
    }
    modes_list <- .getResources(
        eh, modes_metadat[fileIdx, c("Title", "RDataPath")], prefix, verbose
    )
    names(modes_list) <- gsub(prefix, "", names(modes_list))

    eh_experiments <- ExperimentList(modes_list)[resultModes]

    ess_names <- c("colData", "metadata", "sampleMap")

    ess_idx <- .conditionToIndex(ess_names, eh_assays,
        function(x) grepl(x, eh_assays))

    ess_list <- .getResources(eh,
        modes_metadat[ess_idx, c("Title", "RDataPath")], prefix, verbose)
    names(ess_list) <- gsub(prefix, "", names(ess_list))

    c(list(experiments = eh_experiments), ess_list)
}

.mergeLowColData <- function(x) {
    newcoldata <- Reduce(
        function(x, y) {
            S4Vectors::merge(x, y, by = "row.names", all = TRUE)
        },
        lapply(x, colData)
    )
    if (length(x) > 1L) {
        rownames(newcoldata) <- newcoldata[["Row.names"]]
        newcoldata <- newcoldata[, -which(colnames(newcoldata) == "Row.names")]
    }
    newcoldata
}
