.removeExt <- function(fnames, ext = "\\.[Rr][Dd][Aa]$") {
    gsub(ext, "", fnames)
}

.modesAvailable <- function(listfiles) {
    slots <- c("metadata", "colData", "sampleMap")
    modes <- gsub("(^[a-z]*_)(.*)", "\\2", listfiles)
    modes <- .removeExt(modes)
    sort(modes[!modes %in% slots])
}

#' Mouse Gastrulation Multi-modal Data
#'
#' @description scNMT assembles data on-the-fly from `ExperimentHub` to
#'     provide a \linkS4class{MultiAssayExperiment} container. The
#'     `mouse_gastrulation` dataset is obtained from
#'     \insertCite{Argelaguet2019-et;textual}{SingleCellMultiModal}
#'
#' @param dataType character(1) Indicates study that produces this type of
#'     data (default: 'mouse_gastrulation')
#'
#' @param modes character() The assay types or modes of data to obtain these
#'     include single cell Chromatin Accessibilty ("acc"), Methylation ("met"),
#'     RNA-seq ("rna") by default.
#'
#' @param dry.run logical(1) Whether to return the dataset names before actual
#'     download (default TRUE)
#'
#' @param verbose logical(1) Whether to show the dataset currently being
#'     (down)loaded (default TRUE)
#'
#' @seealso SingleCellMultiModal-package
#'
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment}
#'
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/}
#'
#' @references
#'     \insertCite{Argelaguet2019-et}{SingleCellMultiModal}
#'
#' @examples
#' scNMT(dataType = "mouse_gastrulation", modes = "*", dry.run = TRUE)
#'
#' @export scNMT
scNMT <-
    function(
        dataType = "mouse_gastrulation", modes = "*", dry.run = TRUE,
        verbose = TRUE, ...
    )
{
    modes_file <- system.file("extdata", "metadata.csv",
        package = "SingleCellMultiModal", mustWork = TRUE)

    dataType <- tolower(dataType)
    stopifnot(is.character(dataType), length(dataType) == 1L, !is.na(dataType))

    modes_metadat <- read.csv(modes_file, stringsAsFactors = FALSE)
    modes_metadat <- modes_metadat[modes_metadat[["dataType"]] == dataType, ]
    eh_assays <- modes_metadat[["ResourceName"]]
    modesAvail <- .modesAvailable(eh_assays)
    if (identical(modes, "*") && dry.run) {
        message("See the list below for available datasets for")
        cat(dataType, ":\n",
            paste(
                strwrap(paste(modesAvail, collapse = " "), width = 46),
            collapse = "\n "),
        "\n")
        return(invisible())
    }
}
