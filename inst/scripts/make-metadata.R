allextpat <- "\\.[Rr][Dd][Aa]$"

.get_Description <- function(data_name, dataType) {
    paste(data_name, "data specific to the", toupper(dataType), "project")
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

setwd("~/github/SingleCellMultiModal")
source("inst/extdata/docuData/scNMT.R")

make_metadata <- function(
    directory = "~/data/scmm",
    dataDirs = "mouse_gastrulation",
    ext_pattern = "\\.[Rr][Dd][Aa]$",
    doc_file = "inst/extdata/docuData/singlecellmultimodal.csv",
    pkg_name = "SingleCellMultiModal",
    dry.run = TRUE,
    append = FALSE)
{
    if (!identical(basename(getwd()), pkg_name))
        stop("Run 'make_metadata()' from directory: ", pkg_name)

    exdata <- "inst/extdata"

    if (!dir.exists(exdata))
        dir.create(exdata)

    metafile <- file.path(exdata, "metadata.csv")

    metadat <- MetaHubCreate(
        base_dir = directory,
        data_dirs = dataDirs,
        ext_pattern = ext_pattern,
        doc_file = doc_file,
        pkg_name = pkg_name
    )

    if (!dry.run) {
        file.remove(metafile)
        readr::write_csv(metadat, metafile, append = append, col_names = TRUE)
    }

    metadat
}

