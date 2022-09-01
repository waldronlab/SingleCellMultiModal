############################################################
#
# author: Ludwig Geistlinger
# date: 2021-03-24 18:17:27
#
# descr: G&T-seq data retrieval
#
############################################################

#' Parallel sequencing data of single-cell genomes and transcriptomes
#'
#' @description GTseq assembles data on-the-fly from `ExperimentHub` to
#'     provide a \linkS4class{MultiAssayExperiment} container. The `DataType`
#'     argument provides access to the `mouse_embryo_8_cell` dataset as obtained
#'     from Macaulay et al. (2015). Protocol information for this dataset is
#'     available from Macaulay et al. (2016). See references.
#'
#' @details G&T-seq is a combination of Picoplex amplified gDNA sequencing
#'   (genome) and SMARTSeq2 amplified cDNA sequencing (transcriptome) of the
#'   same cell. For more information, see Macaulay et al. (2015).
#'   \itemize{
#'       \item{mouse_embryo_8_cell:}
#'       \itemize{
#'           \item{genomic} - integer copy numbers as detected from scDNA-seq
#'           \item{transcriptomic} - raw read counts as quantified from scRNA-seq
#'       }
#'   }
#'
#' @section metadata:
#'   The `MultiAssayExperiment` metadata includes the original function call
#'   that saves the function call and the data version requested.
#'
#' @param DataType character(1) Indicates study that produces this type of
#'   data (default: 'mouse_embryo_8_cell')
#'
#' @param modes character() A wildcard / glob pattern of modes, such as
#'   \code{"*omic"}. A wildcard of \code{"*"} will return all modes including
#'   copy numbers ("genomic") and RNA-seq read counts ("transcriptomic"),
#'   which is the default.
#'
#' @param version character(1). Currently, only version '1.0.0'.
#'
#' @param dry.run logical(1) Whether to return the dataset names before actual
#'   download (default TRUE)
#'
#' @param verbose logical(1) Whether to show the dataset currently being
#'   (down)loaded (default TRUE)
#'
#' @param ... Additional arguments passed on to the
#'   \link[ExperimentHub]{ExperimentHub-class} constructor
#'
#' @seealso SingleCellMultiModal-package
#'
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment} or
#'     informative `data.frame` when `dry.run` is `TRUE`
#'
#' @source \url{https://www.ebi.ac.uk/ena/browser/view/PRJEB9051}
#'
#' @references
#'   Macaulay et al. (2015) G&T-seq: parallel sequencing of single-cell
#'   genomes and transcriptomes. Nat Methods, 12:519–22.
#'
#'   Macaulay et al. (2016) Separation and parallel sequencing of the genomes
#'   and transcriptomes of single cells using G&T-seq. Nat Protoc, 11:2081–103.
#'
#' @md
#'
#' @examples
#'
#' GTseq()
#'
#' @export GTseq
GTseq <-
    function(
        DataType = "mouse_embryo_8_cell", modes = "*",
        version = "1.0.0", dry.run = TRUE, verbose = TRUE, ...
    )
{
    stopifnot(.isSingleChar(version), .isSingleChar(DataType))
    meta <- list(call = match.call())

    ess_list <- .getResourcesList(
        prefix = "GTseq_",
        datatype = DataType,
        modes = modes,
        version = version,
        dry.run = dry.run,
        verbose = verbose,
        ...
    )

    if (dry.run) { return(ess_list) }

    cdat <- ess_list[["colData"]]
    prim.ids <- rep(paste0("cell", seq_len(112)), 2)
    smap <- S4Vectors::DataFrame(
        assay = tolower(cdat[,"Comment.LIBRARY_SOURCE."]),
        primary = prim.ids,
        colname = cdat[,"Sample.ID"]
    )

    rcols <- c("organism", "sex", "cell.type")
    rcols <- paste0("Characteristics.", rcols, ".")
    cdat <- cdat[seq_len(112), rcols]
    rownames(cdat) <- prim.ids[seq_len(112)]

    MultiAssayExperiment(
        experiments = ess_list[["experiments"]],
        colData = cdat,
        sampleMap = smap,
        metadata = c(meta, as.list(ess_list[["metadata"]]))
    )
}
