setwd("~/github/SingleCellMultiModal")
source("inst/scripts/make-metadata.R")

.inferSource <- function(filepaths) {
    lfiles <- strsplit(filepaths, "\\.")
    exts <- mapply(`[`, lfiles, lengths(lfiles))
    uexts <- toupper(exts)
    vTypes <- AnnotationHubData::getValidSourceTypes()
    uTypes <- toupper(vTypes)
    allvalid <- all(uexts %in% uTypes)
    if (!allvalid)
        stop("Source types not supported: ", paste0(exts[!allvalid],
            collapse = ", "), "\n See 'AnnotationHubData::getValidSources()'",
            call. = FALSE)
    vTypes[match(uexts, uTypes)]
}

.stdLength <- function(metalist, replength) {
    lapply(metalist, function(field) {
        if (length(field) == 1L)
            rep(field, replength)
        else
            field
    })
}

.loadRDAList <- function(filepaths, ext_pat) {
    basefiles <- gsub(ext_pat, "", basename(filepaths))
    filepaths <- setNames(filepaths, basefiles)
    Map(function(filepath, name) {
        OBJENV <- new.env(parent = emptyenv())
        load(filepath, envir = OBJENV)
        OBJENV[[name]]
    }, filepaths, basefiles)
}


## alist() with formals()<-
## fancyFUN <- function() {}
## formals(fancyFUN) <- alist()

MetaHubCreate <-
    function(base_dir, data_dirs, ext_pattern, doc_file, pkg_name)
{
    locations <- file.path(base_dir, data_dirs)
    stopifnot(
        dir.exists(base_dir), all(dir.exists(locations)),
        is.character(ext_pattern), !is.na(ext_pattern),
        identical(length(ext_pattern), 1L),
        file.exists(doc_file), is.character(doc_file), !is.na(doc_file),
        identical(length(doc_file), 1L)
    )
    fpathlist <- lapply(locations, function(locs) {
        list.files(
            locs, pattern = ext_pattern, full.names = TRUE, recursive = TRUE
        )
    })
    docFrame <- read.table(doc_file, header = TRUE)
    docList <- split(docFrame, seq_len(nrow(docFrame)))
    dataTypes <- data_dirs
    replengths <- lengths(fpathlist)
    namelist <- lapply(fpathlist, basename)

    metaList <- Map(
        function(dataType, doc_file, resnames, filepaths, replength) {
            message("Working on: ", basename(dataType))
            dataList <- .loadRDAList(filepaths, ext_pattern)
            hubmeta <- R6::R6Class("EHubMeta",
                public = list(
                    Title = NA_character_,
                    Description = NA_character_,
                    BiocVersion = as.character(BiocManager::version()),
                    Genome = character(1L),
                    SourceType = NA_character_,
                    SourceUrl = character(1L),
                    SourceVersion = NA_character_,
                    Species = character(1L),
                    TaxonomyId = character(1L),
                    Coordinate_1_based = NA,
                    DataProvider = character(1L),
                    Maintainer = NA_character_,
                    RDataClass = NA_character_,
                    DispatchClass = .get_DispatchClass(resnames, ext_pattern),
                    Location_Prefix = NA_character_,
                    RDataPath = NA_character_,
                    ResourceName = resnames,
                    DataType = dataType,

                    initialize = function(doc_file)
                    {
                        if (is.na(self$Title))
                            self$Title <- gsub(ext_pattern, "",
                                basename(filepaths))
                        if (is.na(self$Description))
                            self$Description <- paste(self$Title,
                                "data specific to the", toupper(self$DataType),
                                "project")
                        if (is.na(self$SourceType))
                            self$SourceType <- .inferSource(filepaths)
                        if (is.na(self$SourceVersion))
                            self$SourceVersion <- "1.0.0"
                        if (is.na(self$Maintainer))
                            self$Maintainer <- utils::maintainer(pkg_name)
                        if (is.na(self$RDataClass))
                            self$RDataClass <- .getRDataClass(dataList)
                        if (is.na(self$Location_Prefix))
                            self$Location_Prefix <- NULL
                        if (is.na(self$RDataPath))
                            self$RDataPath <- file.path(pkg_name, self$DataType,
                                self$ResourceName)
                        lapply(names(doc_file), function(i) {
                            assign(i, doc_file[[i]], self)
                        })
                    },
                    generate = function() {
                        lnames <- !names(self) %in%
                            c(".__enclos_env__", "clone", "generate",
                                "initialize")
                        initList <- mget(names(self)[lnames], envir = self)
                        initList <- Filter(function(x) !is.null(x), initList)
                        flist <- .stdLength(initList, replength)
                        do.call(data.frame, c(flist, stringsAsFactors = FALSE))
                    }
                )
            )
            nhub <- hubmeta$new(doc_file)
            nhub$generate()
    }, dataType = dataTypes, doc_file = docList, resnames = namelist,
    filepaths = fpathlist, replength = replengths)
    do.call(
        function(...) {
            rbind.data.frame(..., make.row.names = FALSE,
                stringsAsFactors = FALSE)
        },
    metaList)
}
