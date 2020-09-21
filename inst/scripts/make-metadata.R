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

make_metadata <- function(
    directory = "~/data/scmm",
    dataDirs = c(rep("mouse_gastrulation", 2), "mouse_visual_cortex"),
    version = c("1.0.0", "2.0.0"),
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
