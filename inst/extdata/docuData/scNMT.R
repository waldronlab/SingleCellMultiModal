# TODO: create file and use within constructor function
scnmtData <- data.frame(
    DataProvider = "Dept. of Bioinformatics, The Babraham Institute, United Kingdom",
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = "https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ",
    SourceType = "RDS",
    stringsAsFactors = FALSE
)

## read static file as input to functions

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

library(R6)

## alist() with formals()<-
## fancyFUN <- function() {}
## formals(fancyFUN) <- alist()

MetaHubCreate <- function(base_dir, data_dir, ext_pattern, doc_file, data_list,
    pkg_name = "SingleCellMultiModal")
{
    stopifnot(
        dir.exists(base_dir), dir.exists(data_dir),
        is.character(ext_pattern), !is.na(ext_pattern),
        identical(length(ext_pattern), 1L),
        file.exists(doc_file), is.character(doc_file), !is.na(doc_file),
        identical(length(doc_file), 1L)
    )
    location <- file.path(base_dir, data_dir)
    filepaths <- list.files(
        location, pattern = ext_pattern, full.names = TRUE, recursive = TRUE
    )
    dataType <- data_dir
    replength <- length(filepaths)
    resnames <- basename(filepaths)

    R6::R6Class("EHubMeta",
        public = list(
            Title = NA_character_,
            Description = NA_character_,
            BiocVersion = as.character(BiocManager::version()),
            Genome = character(1L),
            SourceType = NA_character_,
            SourceUrl = character(1L),
            SourceVersion = "1.0.0",
            Species = character(1L),
            TaxonomyId = character(1L),
            Coordinate_1_based = NA,
            DataProvider = character(1L),
            Maintainer = NA_character_,
            RDataClass = NA_character_,
            DispatchClass = .get_DispatchClass(resnames),
            RDataPath = NA_character_,
            ResourceName = resnames,
            DataType = dataType,

            initialize = function(Title, Description, BiocVersion, Genome,
                SourceType, SourceUrl, SourceVersion, Species, TaxonomyId,
                Coordinate_1_based, DataProvider, Maintainer, RDataClass,
                DispatchClass, ResourceName, RDataPath, DataType)
            {
                if (is.na(Title))
                    self$Title <- gsub(ext_pattern, "", basename(filepaths))
                if (is.na(Description))
                    self$Description <- paste(Title, "data specific to the",
                        toupper(DataType), "project")
                if (is.na(SourceType))
                    self$Description <- .inferSource(filepaths)
                if (is.na(Maintainer))
                    self$Maintainer <- utils::maintainer(pkg_name)
                if (is.na(RDataClass))
                    self$RDataClass <- .getRDataClass(dataList)
                if (is.na(RDataPath))
                    self$RDataPath <- file.path(pkg_name, dataType, ResourceName)
            },
            generate = function() {
                # spill guts to data.frame

            }
        )
    )
}
