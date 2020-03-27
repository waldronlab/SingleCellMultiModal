# upload files to AWS S3
allextpat <- "\\.[Rr][Dd][Aa]$"

upload_aws <- function(
    dataList, dataType = "mouse_gastrulation", directory = "~/data/scmm",
    upload = FALSE, fileExt = allextpat
) {
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
