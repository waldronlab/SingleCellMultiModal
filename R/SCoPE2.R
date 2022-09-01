#' Single-cell RNA sequencing and proteomics
#'
#' @description SCoPE2 assembles data on-the-fly from `ExperimentHub`
#'   to provide a \linkS4class{MultiAssayExperiment} container. The
#'   `DataType` argument provides access to the `SCoPE2` dataset as
#'   provided by Specht et al. (2020; DOI: http://dx.doi.org/10.1101/665307).
#'   The article provides more information about the data
#'   acquisition and pre-processing.
#'
#' @details The SCoPE2 study combined scRNA-seq (transcriptome) and
#'   single-cell proteomics. The cells are monocytes that undergo
#'   macrophage differentiation. No annotation is available for the
#'   transcriptome data, but batch and cell type annotations are
#'   available for the proteomics data. The transcriptomics and
#'   proteomics data were not measured from the same cells but from
#'   a distinct set of cell cultures.
#'
#'   \itemize{
#'       \item{SCoPE2:}
#'       \itemize{
#'           \item{scRNAseq1} - single-cell transcriptome (batch 1)
#'           \item{scRNAseq2} - single-cell transcriptome (batch 2)
#'           \item{scp} - single-cell proteomics
#'       }
#'   }
#'
#' @inheritParams scNMT
#'
#' @param DataType character(1) Indicates study that produces this type of
#'   data (default: 'macrophage_differentiation')
#'
#' @param modes character() A wildcard / glob pattern of modes, such as
#'   `"rna"`. A wildcard of `"*"` will return all modes, that are
#'   transcriptome ("rna") or proteome ("protein") which is the
#'   default.
#'
#' @param version character(1), currently only version '1.0.0' is
#'   available
#'
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment} or
#'   informative `data.frame` when `dry.run` is `TRUE`
#'
#' @seealso SingleCellMultiModal-package
#'
#' @source All files are linked from the slavovlab website
#'     \url{https://scope2.slavovlab.net/docs/data}
#'
#' @references
#'   Specht, Harrison, Edward Emmott, Aleksandra A. Petelski, R.
#'   Gray Huffman, David H. Perlman, Marco Serra, Peter Kharchenko,
#'   Antonius Koller, and Nikolai Slavov. 2020. “Single-Cell
#'   Proteomic and Transcriptomic Analysis of Macrophage
#'   Heterogeneity.” bioRxiv. https://doi.org/10.1101/665307.
#'
#' @md
#'
#' @examples
#'
#' SCoPE2(DataType = "macrophage_differentiation",
#'        modes = "*",
#'        version = "1.0.0",
#'        dry.run = TRUE)
#'
#' @export
SCoPE2 <- function(
    DataType = "macrophage_differentiation",
    modes = "*",
    version = "1.0.0",
    dry.run = TRUE,
    verbose = TRUE,
    ...
) {
    if (version != "1.0.0")
        stop("Only version '1.0.0' is available.")

    ## Retrieve the different resources from ExperimentHub
    ess_list <- .getResourcesList(
        prefix = "macrophage_",
        datatype = DataType,
        modes = modes,
        version = version,
        dry.run = dry.run,
        verbose = verbose,
        ...
    )
    ## If dry.run, return only the information table
    if (dry.run) return(ess_list)
    ## Get the colData
    cd <- .mergeLowColData(ess_list[["experiments"]])
    colnames(cd)[which(colnames(cd) == "Batch")] <- "batch_Chromium"

    ## Construct and return the MAE object
    MultiAssayExperiment(
        experiments = ess_list[["experiments"]],
        colData = cd
    )
}
