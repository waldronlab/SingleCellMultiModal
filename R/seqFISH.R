#' Single-cell spatial + Gene Expression
#'
#' @description seqFISH function assembles data on-the-fly from `ExperimentHub`
#'     to provide a \linkS4class{MultiAssayExperiment} container. Actually
#'     the `DataType` argument provides access to the available datasets
#'     associated to the package.
#' @details seq FISH data are a combination of single cell spatial coordinates
#'     and transcriptomics for a few hundreds of genes.
#'     seq-FISH data can be combined for example with scRNA-seq data to unveil
#'     multiple aspects of cellular behaviour based on their spatial
#'     organization and transcription.
#'
#'     Available datasets are:
#'     \itemize{
#'         \item{mouse_visual_cortex: } combination of seq-FISH data as obtained
#'         from Zhu et al. (2018) and scRNA-seq data as obtained from
#'         Tasic et al. (2016),
#'         Version 1.0.0 returns the full scRNA-seq data matrix, while version
#'         2.0.0 returns the processed and subsetted scRNA-seq data matrix
#'         (produced for the Mathematical Frameworks for Integrative Analysis
#'         of Emerging Biological Data Types 2020 Workshop)
#'         The returned seqFISH data are always the processed ones for the same
#'         workshop.
#'         Additionally, cell types annotations are available in the `colData`
#'         through the `class` column in the seqFISH `assay`.
#'         \itemize{
#'             \item{scRNA_Counts} - Tasic scRNA-seq gene count matrix
#'             \item{scRNA_Labels} - Tasic scRNA-seq cell labels
#'             \item{seqFISH_Coordinates} - Zhu seq-FISH spatial coordinates
#'             \item{seqFISH_Counts} - Zhu seq-FISH gene counts matrix
#'             \item{seqFISH_Labels} - Zhu seq-FISH cell labels
#'         }
#'     }
#'
#' @inheritParams scNMT
#'
#' @param DataType character(1) indicating the identifier of the dataset to
#'     retrieve.  (default "mouse_visual_cortex")
#'
#' @param modes character( ) The assay types or modes of data to obtain these
#'     include seq-FISH and scRNA-seq data by default.
#'
#' @return A \linkS4class{MultiAssayExperiment} of seq-FISH data
#'
#' @author Dario Righelli <dario.righelli <at> gmail.com>
#'
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#'
#' @examples
#'
#' seqFISH(DataType = "mouse_visual_cortex", modes = "*", version = "2.0.0",
#'     dry.run = TRUE)
#'
#' @export
seqFISH <-
    function(
        DataType="mouse_visual_cortex", modes="*", version,
        dry.run=TRUE, verbose=TRUE, ...
    )
{
    ess_list <- .getResourcesList(prefix = "seqfish_", datatype = DataType,
        modes = modes, version = version, dry.run = dry.run,
        verbose = verbose, ...)

    if (dry.run) { return(ess_list) }

    modes_list <- ess_list[["experiments"]]

    switch(DataType,
        "mouse_visual_cortex" = {
            mae <- .mouse_visual_cortex(modes_list=modes_list,
                                        version=version)
        },
        ## Add here other seqFISH datasets based on DataType identifier
        {
            stop("Unrecognized seqFISH dataset name")
        }
    )

    return(mae)
}

.mouse_visual_cortex <- function(modes_list, version)
{
    res <- paste0("scRNA",
        if (identical(version, "1.0.0")) "_Full" else "",
        "_", c("Counts", "Labels")
    )

    ## discrepancy between labels in counts and colData
    counts <- as.matrix(modes_list[[res[1]]])
    ## rowData is duplicate of rownames [removed]
    coldata <- modes_list[[res[2]]]
    vIDs <- intersect(rownames(coldata), colnames(counts))
    counts <- counts[, vIDs]
    coldata <- coldata[vIDs, ]

    sce <- SingleCellExperiment::SingleCellExperiment(
        colData=coldata,
        assays=S4Vectors::SimpleList(counts=counts)
    )

    se <- SpatialExperiment::SpatialExperiment(
        rowData=rownames(modes_list$seqFISH_Counts),
        colData=modes_list$seqFISH_Labels,
        assays=S4Vectors::SimpleList(
            counts=as.matrix(modes_list$seqFISH_Counts)),
        spatialData=DataFrame(modes_list$seqFISH_Coordinates),
        spatialCoordsNames=c("x", "y"))

    MultiAssayExperiment(
        experiments = list(seqFISH = se, scRNAseq = sce)
    )
}
