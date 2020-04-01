.checkSize <- function(access, directory = ".", verbose = TRUE) {
    finfo <- getGEOSuppFiles(access, fetch_files = FALSE)
    fur <- as.character(finfo[['url']])

    header <- httr::HEAD(fur)$headers
    header_bytes <- as.numeric(header$`content-length`)

    fpath <- file.path(directory, access, finfo[['fname']])
    local_bytes <- file.size(fpath)

    if (verbose)
        message("url: ", header_bytes, " vs. local: ", local_bytes)
    identical(header_bytes, local_bytes)
}

.getAllSeries <- function(gseacc) {
    if (!is(gseacc, "ExpressionSet"))
        stop("Provide an 'ExpressionSet' for obtaining series from 'colnames'")
    cells <- colnames(gseacc)
    metalist <- lapply(cells, function(x)
        Meta(getGEO(x, GSEMatrix = FALSE))
    )
    allseries <- unique(unlist(
        lapply(metalist, `[[`, "series_id")
    ))
}

.downloadSeries <- function(series, directory = ".") {
    directory <- normalizePath(directory)
    vapply(series, function(x) {
        finfo <- getGEOSuppFiles(x, fetch_files = FALSE)
        fname <- finfo[['fname']]
        fpath <- file.path(directory, x, fname)
        if (file.exists(fpath) && .checkSize(x, directory, verbose = FALSE))
            message("File exists: ", basename(fpath))
        else
            getGEOSuppFiles(x, baseDir = directory)
        fpath
    }, character(1L))
}

.nametodframe <- function(fnames, splitter = "_") {
    filen <- fnames
    fnames <- strsplit(fnames, splitter)
    datf <- lapply(fnames, function(fname) {
        fname[length(fname)] <-
            gsub("\\.tsv\\.gz", "", fname[length(fname)])
        samp <- fname[[1L]]
        inst <- fname %in% c("new", "repeat")
        instance <- ifelse(any(inst), fname[inst], NA_character_)
        fname <- fname[!inst]
        plates <- c(1L:(length(fname)-3L), length(fname))
        expro <- c(1L, (length(fname)-2):length(fname))
        explay <- c(1L, length(fname))

        plate <- paste0(fname[-plates], collapse = "_")
        extract <- paste0(fname[-expro], collapse = "_")
        center <- paste0(fname[-explay], collapse = "_")
        c(samp, extract, instance, plate, center)
    })
    sampmap <- do.call(rbind.data.frame, c(datf, stringsAsFactors = FALSE))
    names(sampmap) <- c("geo_accession", "extraction", "instance", "plate",
        "ext_plate")
    sampmap
}

.splitselect <- function(char, indx = 1L, sep = "_") {
    vapply(strsplit(char, "_"), `[[`, character(1L), indx)
}

.selectInRow <- function(dataframe, term, outcol, colname = NULL) {
    if (!is.null(colname))
        unlist(dataframe[unlist(dataframe[[colname]] == term), outcol])
    else
        unlist(dataframe[term, outcol])
}

.loadEnvObj <- function(filepath, name) {
    OBJENV <- new.env(parent = emptyenv())
    load(filepath, envir = OBJENV)
    object <- OBJENV[[name]]
    object
}

.loadRDAList <- function(metadata_frame) {
    objnames <- .selectInRow(metadata_frame,
        metadata_frame[["experimentFiles"]], "objectNames")

    rdafiles <- .selectInRow(metadata_frame,
        metadata_frame[["experimentFiles"]], "files")

    dataList <- lapply(rdafiles, function(dpath) {
        oname <- .selectInRow(metadata_frame, dpath, "objectNames", "files")
        .loadEnvObj(dpath, oname)
    })
    names(dataList) <- objnames
    dataList
}

