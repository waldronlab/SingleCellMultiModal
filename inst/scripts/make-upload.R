.getDataFiles <- function(directory = "~/data/scmm",
    dataDir = "mouse_gastrulation", pattern = allextpat) {
    location <- file.path(directory, dataDir)
    list.files(location, pattern = pattern, full.names = TRUE, recursive = TRUE)
}

# upload files to AWS S3
allextpat <- "\\.[Rr][Dd][Aa]$"

# IMPORTANT!
# Make sure that AWS_DEFAULT_REGION, AWS_ACCESS_KEY_ID, and
# AWS_SECRET_ACCESS_KEY are set in the ~/.Renviron file

source("make-metadata.R")

upload_aws <- function(
    dataType, directory = "~/data/scmm",
    upload = FALSE, fileExt = allextpat
) {
    if (missing(dataType))
        stop("Enter a 'dataType' folder")
    datafilepaths <- .getDataFiles(
        directory = directory, dataDir = dataType, pattern = fileExt
    )
    bucketLocation <-
        file.path("experimenthub", "SingleCellMultiModal", dataType)
    if (upload)
        AnnotationHubData:::upload_to_S3(file = datafilepaths,
            remotename = basename(datafilepaths),
            bucket = bucketLocation)
}

# upload_aws(dataType = "mouse_gastrulation", upload=TRUE)
# upload_aws(dataType = "mouse_visual_cortex", upload=TRUE)
