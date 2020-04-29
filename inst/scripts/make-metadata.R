allextpat <- "\\.[Rr][Dd][Aa]$"

.getDataFiles <-
function(directory = "~/data/scmm",
    dataDir = "mouse_gastrulation", pattern = allextpat) {
    location <- file.path(directory, dataDir)
    list.files(location, pattern = pattern, full.names = TRUE, recursive = TRUE)
}

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

.makeMetaDF <- function(filepaths, includeSlots = FALSE) {
    namespat <- "^[a-z]*_(.*)"

    basefiles <- gsub(allextpat, "", basename(filepaths))
    obj_slots <-
        if (!includeSlots)
            c("metadata", "colData", "sampleMap")
        else
            NULL

    dfr <- S4Vectors::DataFrame(files = as(filepaths, "List"),
        objectNames = basefiles,
        dataNames = gsub(namespat, "\\1", basefiles),
        dataTypes = vapply(
            strsplit(basefiles, "[_-]"), `[[`, character(1L), 2L)
   )
    dfr[["experimentFiles"]] <- !dfr[["dataTypes"]] %in% obj_slots

    dfr
}

.get_DispatchClass <- function(resource_files) {
    ext_map <- data.frame(
        ext_pattern = allextpat,
        Dispatch = "Rda",
        stringsAsFactors = FALSE
    )
    hitMatrix <- vapply(ext_map[["ext_pattern"]],
        function(pat) grepl(pat, resource_files),
            logical(length(resource_files)))
    ext_map[["Dispatch"]][apply(hitMatrix, 1L, which)]
}

make_metadata <- function(
    directory = "~/data/scmm/",
    dataDirs = "mouse_gastrulation",
    ext_pattern = "\\.[Rr][Dd][Aa]$",
    doc_file = "inst/extdata/docuData/singlecellmultimodal.csv",
    pkg_name = "SingleCellMultiModal",
    append = FALSE)
{
    if (!identical(basename(getwd()), pkg_name))
        stop("Run 'make_metadata()' from directory: ", pkg_name)

    exdata <- "inst/extdata"
    metafile <- file.path(exdata, "metadata.csv")

    if (!dir.exists(exdata))
        dir.create(exdata)

    if (file.exists(metafile))
        file.remove(metafile)

    metadat <- MetaHubCreate(
        base_dir = directory,
        data_dirs = dataDirs,
        ext_pattern = ext_pattern,
        doc_file = doc_file,
        pkg_name = pkg_name
    )

    readr::write_csv(metadat, metafile, append = append, col_names = TRUE)
}

