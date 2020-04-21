.removeExt <- function(fnames, ext = "\\.[Rr][Dd][Aa]$") {
    gsub(ext, "", fnames)
}

.modesAvailable <- function(listfiles) {
    slots <- c("metadata", "colData", "sampleMap")
    modes <- gsub("(^[a-z]*_)(.*)", "\\2", listfiles)
    modes <- .removeExt(modes)
    sort(modes[!modes %in% slots])
}

.searchFromInputs <- function(glob, searchFields) {
    regGlob <- glob2rx(unique(glob))
    res <- unlist(lapply(regGlob, function(x) {
        grep(x, searchFields, ignore.case = TRUE, value = TRUE)
        }))
    if (!length(res))
        stop("No matches found, modify search criteria")
    res
}

.conditionToIndex <- function(startVec, testVec, FUN) {
    logmat <- vapply(startVec, FUN, logical(length(testVec)))
    apply(logmat, 1L, any)
}

.getResources <- function(ExperimentHub, resTable, verbose) {
    fileNames <- stats::setNames(resTable[["RDataPath"]], resTable[["Title"]])
    resources <- lapply(fileNames, function(res) {
        if (verbose)
            message("Working on: ", gsub("\\.rda", "", basename(res)))
        query(ExperimentHub, res)[[1L]]
    })

    resources
}

.test_eh <- function(...) {
    tryCatch({
        ExperimentHub(...)
    }, error = function(e) {
        emsg <- conditionMessage(e)
        if (grepl("Timeout", emsg))
            warning("[experimenthub.bioconductor.org] timeout, localHub=TRUE",
                call.=FALSE)
        ExperimentHub(..., localHub = TRUE)
    })
}


#' Mouse Gastrulation Multi-modal Data
#'
#' @description scNMT assembles data on-the-fly from `ExperimentHub` to
#'     provide a \linkS4class{MultiAssayExperiment} container. The
#'     `mouse_gastrulation` dataset is obtained from Argelaguet et al. (2019)
#'
#' @details scNMT is a combination of RNA-seq (transcriptome) and an adaptation
#'     of Nucleosome Occupancy and Methylation sequencing (NOMe-seq, the
#'     methylome and chromatin accessibility) technologies. For more
#'     information, see Reik et al. (2018) DOI: 10.1038/s41467-018-03149-4
#'
#'     \itemize{
#'         \item{mouse_gastrulation:}
#'         \itemize{
#'             \item{rna} - RNA-seq
#'             \item{acc_*} - chromatin accessibility
#'             \item{met_*} - DNA methylation
#'             \itemize{
#'                 \item{cgi} - CpG islands
#'                 \item{CTCF} - footprints of CTCF binding
#'                 \item{DHS} - DNase Hypersensitive Sites
#'                 \item{genebody} - gene bodies
#'                 \item{p300} - p300 binding sites
#'                 \item{promoter} - gene promoters
#'             }
#'         }
#'     }
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
#' @param ... Additional arguments passed on to the
#'     \link[ExperimentHub]{ExperimentHub-class} constructor
#'
#' @seealso SingleCellMultiModal-package
#'
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment}
#'
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/}
#'
#' @references
#'     Argelaguet et al. (2019)
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
    modes_metadat <- modes_metadat[modes_metadat[["DataType"]] == dataType, ]
    eh_assays <- modes_metadat[["ResourceName"]]
    modesAvail <- .modesAvailable(eh_assays)
    if (identical(modes, "*") && dry.run) {
        message("Available data modes for\n",
            "  ", dataType, ":\n",
            paste(
                strwrap(paste(modesAvail, collapse = ", "),
                    width = 46, indent = 4, exdent = 4),
                collapse = "\n"
            )
        )
        return(invisible())
    }

    resultModes <- .searchFromInputs(modes, modesAvail)
    fileIdx <- .conditionToIndex(resultModes, eh_assays, function(x) grepl(x, eh_assays))
    fileMatches <- modes_metadat[fileIdx, c("Title", "DispatchClass")]

    if (dry.run) { return(fileMatches) }
    eh <- .test_eh(...)
    modes_list <- .getResources(
        eh, modes_metadat[fileIdx, c("Title", "RDataPath")], verbose
    )
    names(modes_list) <- gsub("scnmt_", "", names(modes_list))

    eh_experiments <- ExperimentList(modes_list)

    ess_names <- c("colData", "metadata", "sampleMap")

    ess_idx <- .conditionToIndex(ess_names, eh_assays,
        function(x) grepl(x, eh_assays))

    ess_list <- .getResources(eh,
        modes_metadat[ess_idx, c("Title", "RDataPath")], verbose)
    names(ess_list) <- gsub("scnmt_", "", names(ess_list))

    MultiAssayExperiment(
        experiments = eh_experiments,
        colData = ess_list[["colData"]],
        sampleMap = ess_list[["sampleMap"]],
    )
}
