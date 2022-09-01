#' Single-cell Nucleosome, Methylation and Transcription sequencing
#'
#' @description scNMT assembles data on-the-fly from `ExperimentHub` to
#'     provide a \linkS4class{MultiAssayExperiment} container. The `DataType`
#'     argument provides access to the `mouse_gastrulation` dataset as obtained
#'     from Argelaguet et al. (2019; DOI: 10.1038/s41586-019-1825-8).
#'     Pre-processing code can be seen at
#'     \url{https://github.com/rargelaguet/scnmt_gastrulation}. Protocol
#'     information for this dataset is available at Clark et al. (2018). See
#'     the vignette for the full citation.
#'
#' @details scNMT is a combination of RNA-seq (transcriptome) and an adaptation
#'   of Nucleosome Occupancy and Methylation sequencing (NOMe-seq, the
#'   methylome and chromatin accessibility) technologies. For more
#'   information, see Reik et al. (2018) DOI: 10.1038/s41467-018-03149-4
#'
#'   \itemize{
#'       \item{mouse_gastrulation:}
#'       \itemize{
#'           \item{rna} - RNA-seq
#'           \item{acc_*} - chromatin accessibility
#'           \item{met_*} - DNA methylation
#'           \itemize{
#'               \item{cgi} - CpG islands
#'               \item{CTCF} - footprints of CTCF binding
#'               \item{DHS} - DNase Hypersensitive Sites
#'               \item{genebody} - gene bodies
#'               \item{p300} - p300 binding sites
#'               \item{promoter} - gene promoters
#'           }
#'       }
#'   }
#'   Special thanks to Al J Abadi for preparing the published data in time
#'   for the 2020 BIRS Workshop, see the link here:
#'   url{https://github.com/BIRSBiointegration/Hackathon/tree/master/scNMT-seq}
#'
#' @section versions:
#'   Version '1.0.0' of the scNMT mouse_gastrulation dataset includes all of
#'   the above mentioned assay technologies with filtering of cells based on
#'   quality control metrics. Version '2.0.0' contains all of the cells
#'   without the QC filter and does not contain CTCF binding footprints or
#'   p300 binding sites.
#'
#' @section metadata:
#'   The `MultiAssayExperiment` metadata includes the original function call
#'   that saves the function call and the data version requested.
#'
#' @param DataType character(1) Indicates study that produces this type of
#'   data (default: 'mouse_gastrulation')
#'
#' @param modes character() A wildcard / glob pattern of modes, such as
#'   \code{"acc*"}. A wildcard of \code{"*"} will return all modes including
#'   Chromatin Accessibilty ("acc"), Methylation ("met"), RNA-seq ("rna")
#'   which is the default.
#'
#' @param version character(1) Either version '1.0.0' or '2.0.0' depending on
#'   data version required (default '1.0.0'). See version section.
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
#'   informative `data.frame` when `dry.run` is `TRUE`
#'
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/}
#'
#' @references
#'   Argelaguet et al. (2019)
#'
#' @md
#'
#' @examples
#'
#' scNMT(DataType = "mouse_gastrulation", modes = "*",
#'     version = "1.0.0", dry.run = TRUE)
#'
#' @export scNMT
scNMT <-
    function(
        DataType = "mouse_gastrulation", modes = "*", version = "1.0.0",
        dry.run = TRUE, verbose = TRUE, ...
    )
{
    stopifnot(.isSingleChar(version), .isSingleChar(DataType))
    meta <- list(call = match.call(), version = version)

    if (missing(version) || !version %in% c("1.0.0", "2.0.0"))
        stop("Enter version '1.0.0' or '2.0.0'; see '?scNMT' for details.")

    ess_list <- .getResourcesList(prefix = "scnmt_", datatype = DataType,
        modes = modes, version = version, dry.run = dry.run,
        verbose = verbose, ...)

    if (dry.run) { return(ess_list) }

    MultiAssayExperiment(
        experiments = ess_list[["experiments"]],
        colData = ess_list[["colData"]],
        sampleMap = ess_list[["sampleMap"]],
        metadata = meta
    )
}
