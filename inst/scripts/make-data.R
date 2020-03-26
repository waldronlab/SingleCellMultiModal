# get data from cloudstor
# https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
## ./output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
if (FALSE) {
    library(MultiAssayExperiment)

    mae <- readRDS("scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds")
    ddir <- "~/data/scmm/mouse_gastrulation"

    if (!dir.exists(ddir))
        dir.create(ddir, recursive = TRUE)

    exportClass(mae, ddir, fmt = "csv")
}

# convert .csv files to .rda matrices
.convertData <- function(
    directory = "~/data/scmm/",
    dataDir = "mouse_gastrulation",
    pattern = ".csv")
{
    location <- file.path(directory, dataDir)
    csvs <- list.files(location, pattern = pattern, full.names = TRUE,
        recursive = TRUE)
    invisible(
        lapply(csvs, function(csvfile) {
            objname <- gsub(pattern, "", basename(csvfile))
            readin <- as.data.frame(readr::read_csv(csvfile))
            if (!objname %in% c("scnmt_colData", "scnmt_sampleMap")) {
                readin <- data.matrix(readin)
            }
            rownames(readin) <- readin[, 1L]
            readin <- readin[, -1L]
            rdafile <- gsub("csv", "rda", csvfile)
            assign(objname, readin)
            save(list = objname, file = rdafile)
        })
    )
}

