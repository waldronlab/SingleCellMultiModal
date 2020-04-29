# upload files to AWS S3
allextpat <- "\\.[Rr][Dd][Aa]$"

# IMPORTANT! First do:
# aws sts get-session-token
# then add key value pairs to ~/.Renviron
ab <- system2("aws", c("sts", "get-session-token"), stdout = TRUE)
creds <- strsplit(ab, "\t")[[1]][c(2,4,5)]
credlines <- paste0(
    c("AWS_ACCESS_KEY_ID=", "AWS_SECRET_ACCESS_KEY=", "AWS_SESSION_TOKEN="),
    creds
)
renv <- readLines("~/.Renviron")
renv <- renv[-grep("AWS_[AS]", renv)]
writeLines(renv, con = "~/.Renviron")
write(credlines, file = "~/.Renviron", append = TRUE)

source("make-metadata.R")

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

upload_aws(upload=TRUE)
