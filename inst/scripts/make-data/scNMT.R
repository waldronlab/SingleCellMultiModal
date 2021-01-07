# get data from cloudstor
# https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ/download?path=%2Foutput&files=scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
## ./output/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds
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

