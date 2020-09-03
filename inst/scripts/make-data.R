# get data from cloudstor
# https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
## ./output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
if (FALSE) {
    library(MultiAssayExperiment)

    ddir <- "~/data/scmm/mouse_gastrulation"

    if (!dir.exists(ddir))
        dir.create(ddir, recursive = TRUE)

#   old
#   "scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds"
    scnmt <- readRDS(
        file.path(ddir, "allcells",
            "scnmtseq_gastrulation_mae_AllCells.rds"
        )
    )

    exportClass(scnmt, ddir, fmt = "csv")
}

# convert .csv files to .rda matrices
.convertData <- function(
    directory = "~/data/scmm/",
    dataDir = "mouse_gastrulation",
    pattern = ".csv")
{
    location <- file.path(directory, dataDir)
    csvs <- list.files(location, pattern = pattern, full.names = TRUE,
        recursive = FALSE)
    invisible(
        lapply(csvs, function(csvfile) {
            objname <- gsub(pattern, "", basename(csvfile))
            readin <- as.data.frame(readr::read_csv(csvfile))
            rnames <- readin[[1L]]

            if (!objname %in% c("scnmt_colData", "scnmt_sampleMap"))
                readin <- data.matrix(readin[, -1])
            else if (identical(objname, "scnmt_colData"))
                names(readin)[1] <- "cellID"
            else
                readin <- readin[, -1]

            if (!objname %in% "scnmt_sampleMap")
                rownames(readin) <- rnames

            assign(objname, readin)
            rdafile <- gsub("csv", "rda", csvfile)
            save(list = objname, file = rdafile)
        })
    )
}
