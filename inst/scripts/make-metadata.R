setwd("~/gh/SingleCellMultiModal")

.getSourceType <- function(filepaths) {
    lfiles <- strsplit(basename(filepaths), "\\.")
    exts <- vapply(lfiles,
        function(x) { paste(x[-1], collapse = ".") }, character(1L))
    uexts <- toupper(exts)
    uexts <- gsub("[Hh]5", "HDF5", uexts)
    uexts <- gsub("[Mm][Tt][Xx]\\.[Gg][Zz]", "MTX", uexts)
    vTypes <- AnnotationHubData::getValidSourceTypes()
    uTypes <- toupper(vTypes)
    allvalid <- all(uexts %in% uTypes)
    if (!allvalid)
        stop("Source types not supported: ", paste0(exts[!allvalid],
            collapse = ", "), "\n See 'AnnotationHubData::getValidSources()'",
            call. = FALSE)
    res <- vTypes[match(uexts, uTypes)]
    ## hot fix before AnnotationHubData 1.21.2
    gsub("MTX", "mtx.gz", res, fixed = TRUE)
}

doc_helper <-
    function(
        DataProvider, TaxonomyId, Species, SourceUrl, SourceType, DataType, ...
    )
{
    args <- list(...)
    saf <- args[["stringsAsFactors"]]
    saf <- if(!is.null(saf)) saf else FALSE

    input_vals <- list(
        DataProvider = DataProvider, TaxonomyId = TaxonomyId,
        Species = Species, SourceUrl = SourceUrl,
        SourceType = SourceType, DataType = DataType
    )
    clens <- lengths(input_vals)
    zlen <- !clens
    if (any(zlen))
        stop(
            "Provide values for: ",
            paste(names(input_vals)[zlen], collapse = ", ")
        )

    nonstd <- !clens %in% c(max(clens), 1L)
    if (any(nonstd))
        stop("Lengths of inputs must either be 1 or the max length")

    input_vals[clens == 1L] <- lapply(input_vals[clens == 1L],
        function(x) {
            rep(x, max(clens))
        })

    as.data.frame(input_vals, stringsAsFactors = saf)
}

.stdLength <- function(metalist, replength) {
    lapply(metalist, function(field) {
        if (length(field) == 1L)
            rep(field, replength)
        else
            field
    })
}

.loadRDS <- function(filepath) {
    readRDS(filepath)
}

.loadRDA <- function(filepath) {
    basefile <- gsub("\\.[Rr][Dd][Aa]", "", basename(filepath))
    OBJENV <- new.env(parent = emptyenv())
    load(filepath, envir = OBJENV)
    OBJENV[[basefile]]
}

.loadH5 <- function(filepath) {
    HDF5Array::HDF5Array(filepath, "assay001")
}

.loadMTX.GZ <- function(filepath) {
    HCAMatrixBrowser:::.read_mtx(filepath)
}

.loadDataList <- function(filepaths) {
    recipelist <- list(
        "\\.[Rr][Dd][Aa]" = .loadRDA,
        "\\.[Rr][Dd][Ss]" = .loadRDS,
        "\\.[Hh]5" = .loadH5,
        "\\.[Mm][Tt][Xx]\\.[Gg][Zz]" = .loadMTX.GZ
    )
    hitMatrix <- vapply(names(recipelist),
        function(pat) grepl(pat, filepaths),
        logical(length(filepaths))
    )
    allrecipes <- recipelist[apply(hitMatrix, 1L, which)]
    Map(function(x, y) { x(y) }, x = allrecipes, y = filepaths)
}

any.na <- function(x) {
    any(is.na(x))
}

.get_Description <- function(data_name, DataType) {
    paste(data_name, "data specific to the", toupper(DataType), "project")
}

.getRDataClass <- function(dataList) {
    vapply(dataList, function(dataName) {
            if (is.matrix(dataName))
                "matrix"
            else
                class(dataName)
    }, character(1L))
}

.getDispatchClass <- function(resource_files, ext_pat) {
    ext_map <- data.frame(
        ext_pattern = paste0(
            c("[Rr][Dd][Aa]", "[Rr][Dd][Ss]", "[Hh]5", "[Mm][Tt][Xx]\\.[Gg][Zz]"),
            "$"
        ),
        ## currently MTX DispatchClass recipe unavailable
        Dispatch = c("Rda", "Rds", "H5File", "FilePath"),
        stringsAsFactors = FALSE
    )
    hitMatrix <- vapply(ext_map[["ext_pattern"]],
        function(pat) grepl(pat, resource_files),
            logical(length(resource_files)))
    ext_map[["Dispatch"]][apply(hitMatrix, 1L, which)]
}

## alist() with formals()<-
## fancyFUN <- function() {}
## formals(fancyFUN) <- alist()

MetaHubCreate <-
    function(base_dir, data_dirs, ext_pattern, doc_file, version, pkg_name)
{
    locations <- file.path(base_dir, data_dirs, paste0("v", version))
    stopifnot(
        dir.exists(base_dir), all(dir.exists(locations)),
        is.character(ext_pattern), !is.na(ext_pattern),
        identical(length(ext_pattern), 1L),
        file.exists(doc_file), is.character(doc_file), !is.na(doc_file),
        identical(length(doc_file), 1L), is.character(version)
    )
    fpathlist <- lapply(locations, function(locs) {
        list.files(
            locs, pattern = ext_pattern, full.names = TRUE, recursive = TRUE
        )
    })
    docFrame <- read.csv(doc_file, header = TRUE)
    docList <- split(docFrame,
        list(docFrame[["DataType"]], docFrame[["SourceVersion"]]))
    versions <- version
    DataTypes <- data_dirs
    replengths <- lengths(fpathlist)
    namelist <- lapply(fpathlist, basename)

    metaList <- Map(
        function(DataType, doc_file, resnames, filepaths, replength, version) {
            message("Working on: ", basename(DataType), " v", version)
            hubmeta <- R6::R6Class("EHubMeta",
                public = list(
                    Title = NA_character_,
                    Description = NA_character_,
                    BiocVersion = as.character(BiocManager::version()),
                    Genome = NA_character_,
                    SourceType = NA_character_,
                    SourceUrl = character(1L),
                    SourceVersion = version,
                    Species = character(1L),
                    TaxonomyId = character(1L),
                    Coordinate_1_based = NA,
                    DataProvider = character(1L),
                    Maintainer = NA_character_,
                    RDataClass = NA_character_,
                    DispatchClass = .getDispatchClass(resnames, ext_pattern),
                    Location_Prefix = NA_character_,
                    RDataPath = NA_character_,
                    ResourceName = resnames,
                    DataType = DataType,

                    initialize = function(doc_file)
                    {
                        lapply(names(doc_file), function(i) {
                            assign(i, doc_file[[i]], self)
                        })
                        if (is.na(self$Title))
                            self$Title <- gsub(ext_pattern, "",
                                basename(filepaths))
                        if (is.na(self$Description))
                            self$Description <- paste(self$Title,
                                "data specific to the", toupper(self$DataType),
                                "project")
                        if (any.na(self$SourceType))
                            self$SourceType <- .getSourceType(filepaths)
                        if (any.na(self$SourceVersion))
                            self$SourceVersion <- "1.0.0"
                        if (any.na(self$Maintainer))
                            self$Maintainer <- utils::maintainer(pkg_name)
                        if (any.na(self$RDataClass)) {
                            dataList <- .loadDataList(filepaths)
                            self$RDataClass <- .getRDataClass(dataList)
                        }
                        if (is.na(self$Location_Prefix))
                            self$Location_Prefix <- NULL
                        if (is.na(self$RDataPath))
                            self$RDataPath <- file.path(pkg_name,
                                self$DataType, paste0("v", version),
                                self$ResourceName)
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
                ),
                lock_objects = FALSE
            )
            nhub <- hubmeta$new(doc_file)
            nhub$generate()
    }, DataType = DataTypes, doc_file = docList, resnames = namelist,
    filepaths = fpathlist, replength = replengths, version = versions)
    do.call(
        function(...) {
            rbind.data.frame(..., make.row.names = FALSE,
                stringsAsFactors = FALSE)
        },
    metaList)
}

#' Generate the metadata.csv file from a documentation file
#'
#' This function takes a specific folder structure and generates the
#' metadata.csv file for adding to ExperimentHub.
#'
#' @param directory The base folder for _all_ datasets
#'
#' @param dataDirs character() A vector of folder names contained in directory
#'     that corresponds to each project. For multiple versions, repeat the
#'     name of the folder.
#'
#' @param version character() A vector of subfolder versions that is parallel
#'     to `dataDirs` argument, typically `v1.0.0`.
#'
#' @param ext_pattern character(1) A string that matches files within the
#'     above folders to find the data.
#'
#' @param doc_file character(1) A path to the documentation `data.frame` that
#'     tells the function how to fill in the standard columns for data
#'     annotation, for example `DataProvider`, `TaxonomyId`, etc.
#'
#' @param pkg_name character(1) The name of the current package
#'
#' @param dry.run logical(1) Whether to (over)write the `metadata.csv` file or
#'     return as output.
#'
#' @param append logical(1) Whether to append to the current `metadata.csv`
#'     file
#'
#' @return Saves a file under `/inst/extdata/metadata.csv`
#'
#' @examples
#'
#' make_metadata(
#'     directory = "~/data/scmm",
#'     dataDirs = "mouse_gastrulation",
#'     version = c("1.0.0", "2.0.0"),
#'     doc_file = "inst/extdata/docuData/singlecellmultimodalv2.csv",
#'     dry.run = FALSE
#' )
#'
#' make_metadata(
#'     directory = "~/data/scmm",
#'     dataDirs = c(rep("mouse_gastrulation", 2),
#'         rep("mouse_visual_cortex", 2)),
#'     version = rep(c("1.0.0", "2.0.0"), 2),
#'     ext_pattern = "\\.[Rr][Dd][Aa]$",
#'     doc_file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
#'     pkg_name = "SingleCellMultiModal",
#'     dry.run = TRUE,
#' )
#'
#' make_metadata(
#'     directory = "~/data/scmm",
#'     dataDirs = "pbmc",
#'     version = "1.0.0",
#'     ext_pattern = "\\.[Rr][Dd][AaSs]$|\\.[Mm][Tt][Xx]\\.[Gg][Zz]$",
#'     doc_file = "inst/extdata/docuData/singlecellmultimodalv6.csv",
#'     pkg_name = "SingleCellMultiModal",
#'     dry.run = TRUE,
#' )
#'
#' @md
#'
#' @export

make_metadata <- function(
    directory = "~/data/scmm",
    dataDirs = c(rep("mouse_gastrulation", 2), rep("mouse_visual_cortex", 2), "pbmc"),
    version = c(rep(c("1.0.0", "2.0.0"), 2), "1.0.0"),
    ext_pattern = "\\.[Rr][Dd][AaSs]$|\\.[Mm][Tt][Xx]\\.[Gg][Zz]$|\\.[Hh]5$",
    doc_file,
    pkg_name = "SingleCellMultiModal",
    dry.run = TRUE,
    append = FALSE)
{
    if (!identical(basename(getwd()), pkg_name))
        stop("Run 'make_metadata()' from directory: ", pkg_name)

    exdata <- "inst/extdata"

    if (!dir.exists(exdata))
        dir.create(exdata)

    if (missing(doc_file))
        stop("'doc_file' for generating the metadata is missing")

    metafile <- file.path(exdata, "metadata.csv")

    metadat <- MetaHubCreate(
        base_dir = directory,
        data_dirs = dataDirs,
        ext_pattern = ext_pattern,
        doc_file = doc_file,
        version = version,
        pkg_name = pkg_name
    )

    if (!dry.run) {
        if(!append)
        {
            file.remove(metafile)
        }
        readr::write_csv(metadat, metafile, append = append, na="NA")
    }

    metadat
}

# make_metadata(
#     dataDirs = "mouse_gastrulation",
#     version = "1.0.0",
#     doc_file = "inst/extdata/docuData/singlecellmultimodalv1.csv",
#     dry_run = FALSE
# )
#
# make_metadata(
#     directory="CITEseq/",
#     dataDirs = "cord_blood",
#     version = "1.0.0",
#     doc_file = "inst/extdata/docuData/singlecellmultimodalv5.csv",
#     dry.run = FALSE,
#     append=TRUE
# )

# make_metadata(
#     dataDirs = c(rep("mouse_gastrulation", 2), "mouse_visual_cortex"),
#     version = c("1.0.0", "2.0.0", "1.0.0"),
#     doc_file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
#     dry.run = FALSE
# )

# make_metadata(
#     directory = "~/data/scmm",
#     dataDirs = "peripheral_blood",
#     version = "1.0.0",
#     doc_file = "inst/extdata/docuData/singlecellmultimodalv5.csv",
#     dry.run = FALSE,
#     append = TRUE
# )

# make_metadata(
#     directory = "~/data/scmm",
#     dataDirs = "pbmc_10x",
#     version = "1.0.0",
#     doc_file = "inst/extdata/docuData/singlecellmultimodalv6.csv",
#     dry.run = FALSE,
#     append = TRUE
# )

make_metadata(
    directory = "../.localdata/SingleCellMultiModal/",
    dataDirs = "macrophage_differentiation",
    version = "1.0.0",
    doc_file = "inst/extdata/docuData/singlecellmultimodalv7.csv",
    dry.run = FALSE,
    append = TRUE
)


## request to update Maintainer field in older AH resources
# aq <- AnnotationHub::query(eh, "SingleCellMultiModal")
# aq[aq$maintainer == "Marcel Ramos <marcel.ramos@roswellpark.org>" &
#     grepl("v[12]", aq$rdatapath)]

