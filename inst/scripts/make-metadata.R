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

.getMetadata <- function(
    directory, dataDir, ext_pattern, resource_maintainer, resource_biocVersion)
{
    stopifnot(S4Vectors::isSingleString(directory),
        S4Vectors::isSingleString(dataDir))
    ## loop over datasets in each dataDir
    metasets <- lapply(dataDir, function(dataType) {
        datafilepaths <- .getDataFiles(
            directory = directory, dataDir = dataType, pattern = ext_pattern
        )
        message("Working on: ", basename(dataType))
        dfmeta <- .makeMetaDF(datafilepaths, TRUE)
        dataList <- .loadRDAList(dfmeta)
        replen <- length(datafilepaths)

        ResourceName <- basename(datafilepaths)
        Title <- gsub(ext_pattern, "", ResourceName)
        Description <- .get_Description(Title, dataType)
        BiocVersion <- rep(as.character(resource_biocVersion), replen)
        Genome <- rep("", replen)
        SourceType <- rep("RDS", replen)
        SourceUrl <-
            rep("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ",
                replen)
        SourceVersion <- rep("1.0.0", replen)
        Species <- rep("Mus musculus", replen)
        TaxonomyId <- rep("10090", replen)
        Coordinate_1_based <- rep(as.logical(NA), replen)
        DataProvider <-
            rep("Dept. of Bioinformatics, The Babraham Institute, United Kingdom", replen)
        Maintainer <- rep(resource_maintainer, replen)
        RDataPath <- file.path("singlecellmultimodal", dataType, ResourceName)
        RDataClass <- .getRDataClass(dataList)
        DispatchClass <- .get_DispatchClass(ResourceName)
        DataType <- rep(dataType, replen)
        data.frame(Title, Description, BiocVersion, Genome, SourceType, SourceUrl,
                   SourceVersion, Species, TaxonomyId, Coordinate_1_based,
                   DataProvider, Maintainer, RDataClass, DispatchClass,
                   ResourceName, RDataPath, DataType, stringsAsFactors = FALSE)
    })
    do.call(rbind, metasets)
}

make_metadata <- function(
    directory = "~/data/scmm/",
    dataDir = "mouse_gastrulation",
    ext_pattern = "\\.[Rr][Dd][Aa]$",
    resource_maintainer = utils::maintainer("SingleCellMultiModal"),
    resource_biocVersion = BiocManager::version())
{
    if (!identical(basename(getwd()), "SingleCellMultiModal"))
        stop("Run 'make_metadata()' from directory: 'SingleCellMultiModal'")

    exdata <- "inst/extdata"
    metafile <- file.path(exdata, "metadata.csv")

    if (!dir.exists(exdata))
        dir.create(exdata)

    if (file.exists(metafile))
        file.remove(metafile)

    metadat <- .getMetadata(directory = directory, dataDir = dataDir,
        ext_pattern = ext_pattern, resource_maintainer = resource_maintainer,
        resource_biocVersion = resource_biocVersion)

    readr::write_csv(metadat, "inst/extdata/metadata.csv", col_names = TRUE)
}

