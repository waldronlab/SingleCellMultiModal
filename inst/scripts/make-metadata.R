allextpat <- "\\.[Rr][Dd][Aa]$"

.getDataFiles <-
function(directory = "~/data/scmm/",
    dataDir = "scnmt_gastrulation", pattern = allextpat) {
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

## TODO: Update after here
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
    dataTypeFolders <- dir(file.path(directory, dataDir))
    stopifnot(S4Vectors::isSingleString(directory),
        S4Vectors::isSingleString(dataDir))

    metasets <- lapply(dataTypeFolders, function(dataType) {
        message("Working on: ", dataType)
        datafilepaths <- .getDataFiles(directory = directory,
            dataDir = dataDir, dataTypeFolder = dataType, pattern = ext_pattern)
        dfmeta <- .makeMetaDF(datafilepaths, TRUE)
        dataList <- .loadRDAList(dfmeta)
        dataList <- .addMethylation(dfmeta, dataList)
        replen <- length(datafilepaths)

        ResourceName <- basename(datafilepaths)
        Title <- gsub(ext_pattern, "", ResourceName)
        Description <- .get_Description(Title, dataType)
        BiocVersion <- rep(as.character(resource_biocVersion), replen)
        Genome <- rep("", replen)
        SourceType <- rep("TXT", replen)
        SourceUrl <- rep("http://gdac.broadinstitute.org/", replen)
        SourceVersion <- rep("1.1.38", replen)
        Species <- rep("Homo sapiens", replen)
        TaxonomyId <- rep("9606", replen)
        Coordinate_1_based <- rep(as.logical(NA), replen)
        DataProvider <-
            rep("Eli and Edythe L. Broad Institute of Harvard and MIT", replen)
        Maintainer <- rep(resource_maintainer, replen)
        RDataPath <- file.path("curatedTCGAData", ResourceName)
        RDataClass <- .getRDataClass(dataList)
        DispatchClass <- .get_DispatchClass(ResourceName)
        data.frame(Title, Description, BiocVersion, Genome, SourceType, SourceUrl,
                   SourceVersion, Species, TaxonomyId, Coordinate_1_based,
                   DataProvider, Maintainer, RDataClass, DispatchClass,
                   ResourceName, RDataPath, stringsAsFactors = FALSE)
    })
    do.call(rbind, metasets)
}


make_metadata <- function(
    directory = "~/data/scmm/",
    dataDir = "scnmt_gastrulation",
    ext_pattern = "\\.[Rr][Dd][Aa]$",
    resource_maintainer = utils::maintainer("SingleCellMultiModal"),
    resource_biocVersion = BiocManager::version())
{
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

