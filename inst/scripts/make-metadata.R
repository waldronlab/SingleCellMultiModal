allextpat <- "\\.[Rr][Dd][Aa]$"

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

.get_DispatchClass <- function(resource_files, ext_pat) {
    ext_map <- data.frame(
        ext_pattern = ext_pat,
        Dispatch = "Rda",
        stringsAsFactors = FALSE
    )
    hitMatrix <- vapply(ext_map[["ext_pattern"]],
        function(pat) grepl(pat, resource_files),
            logical(length(resource_files)))
    ext_map[["Dispatch"]][apply(hitMatrix, 1L, which)]
}

# setwd("~/gh/SingleCellMultiModal")
source("inst/extdata/docuData/scNMT.R")

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
#'     directory = "~/data/scmm"
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
#' @md
#'
#' @export
make_metadata <- function(
    directory = "~/data/scmm",
    dataDirs = c(rep("mouse_gastrulation", 2), rep("mouse_visual_cortex", 2)),
    version = rep(c("1.0.0", "2.0.0"), 2),
    ext_pattern = "\\.[Rr][Dd][Aa]$",
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
#     dry_run = FALSE
# )
