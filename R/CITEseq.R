.cord_blood <- function(ess_list)
{
    names(ess_list$experiments) <- gsub("_Counts", "", names(ess_list$experiments))
    mse <- MultiAssayExperiment::MultiAssayExperiment(experiments=(ess_list$experiments))
    return(mse)
}


#' CITEseq
#' @description function assembles data on-the-fly from `ExperimentHub`
#'     to provide a \linkS4class{MultiAssayExperiment} container. Actually
#'     the `dataType` argument provides access to the available datasets
#'     associated to the package.
#' @details CITEseq data are a combination of single cell transcriptomics and
#'     about a hundread of cell surface proteins.
#'
#'     Available datasets are:
#'     \itemize{
#'         \item{cord_blood: } a dataset of single cells of cord blood as
#'         provided in Stoeckius et al. (2017).
#'          \itemize{
#'             \item{scRNA_Counts} - Stoeckius scRNA-seq gene count matrix
#'             \item{scADT} - Stoeckius antibody-derived tags (ADT) data
#'             }
#'      }
#'
#' @param DataType character(1) indicating the identifier of the dataset to
#'     retrieve.  (default "cord_blood")
#'
#' @param modes character( ) The assay types or modes of data to obtain these
#'     include scADT and scRNA-seq data by default.
#' @param version character(1) Either version '1.0.0' depending on
#'     data version required.
#' @param dry.run logical(1) Whether to return the dataset names before actual
#'     download (default TRUE)
#'
#' @param verbose logical(1) Whether to show the dataset currently being
#'     (down)loaded (default TRUE)
#'
#' @param ... Additional arguments passed on to the
#'     \link[ExperimentHub]{ExperimentHub-class} constructor
#'
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment} or
#'     informative `data.frame` when `dry.run` is `TRUE`
#' @references Stoeckius et al. (2017)
#' @export
#'
#' @examples
#'
#' mse <- CITEseq(dry.run=FALSE)
#' experiments(mse)
#'
CITEseq <- function(DataType="cord_blood", modes="*",
                    version="1.0.0", dry.run=TRUE, verbose=TRUE, ...)
{

    ess_list <- .getResourcesList(prefix = "citeseq_", datatype = DataType,
        modes=modes, version=version, dry.run=dry.run, verbose=verbose, ...)
    if (!dry.run) {
        mse <- switch(
            DataType,
            "cord_blood" = { .cord_blood(ess_list=ess_list) },
                ## Add here other CITE-seq datasets based on DataType identifier
            { stop("Unrecognized CITE-seq dataset name") }
        )
        return(mse)
    } else {
        return(ess_list)
    }

}


#' CITEseqMseToSce
#' @description converts a MultiAssayExperiment object with CITEseq data into
#' a SingleCellExperiment object to be used with already known methods and
#' packages in literature.
#'
#' @param mse a MultiAssayExperiment object with scRNA and scADT named
#' experiments
#'
#' @return a SingleCellExperiment object as widely with scRNA data as counts
#' and scADT data as altExps
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom methods is
#' @export
#'
#' @examples
#'
#' mse <- CITEseq(dry.run=FALSE)
#' sce <- CITEseqMseToSce(mse)
#'
CITEseqMseToSce <- function(mse)
{
    stopifnot(is(mse, "MultiAssayExperiment"))

    scrna <- experiments(mse)[[grep("scRNA", names(mse))]]
    scadt <- SummarizedExperiment(experiments(mse)[[grep("scADT", names(mse))]])
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=scrna),
                            altExps=scadt)
    return(sce)
}

