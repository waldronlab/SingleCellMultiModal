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

.queryResources <- function(ExperimentHub, resTable, verbose) {
    fileNames <- stats::setNames(resTable[["RDataPath"]], resTable[["Title"]])
    lapply(fileNames, function(res) {
        if (verbose)
            message("Working on: ", gsub("\\.rda", "", basename(res)))
        query(ExperimentHub, res)
    })
}

.getResources <- function(ExperimentHub, resTable, verbose) {
    infos <- .queryResources(ExperimentHub, resTable, verbose)
    lapply(infos, `[[`, 1L)
}

.getResourceInfo <- function(ExperimentHub, resTable, verbose) {
    infos <- .queryResources(ExperimentHub, resTable, verbose)
    resID <- vapply(infos, names, character(1L))
    restab <- AnnotationHub::getInfoOnIds(ExperimentHub, resID)
    restab <-
        restab[, !names(restab) %in% c("fetch_id", "status", "biocversion")]
    sizes <- as.numeric(restab[["file_size"]])
    class(sizes) <- "object_size"
    restab <- as.data.frame(append(
        restab,
        list(file_size = format(sizes, units = "Mb")),
        which(names(restab) == "title")
    ))
    restab[, -length(restab)]
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

.isSingleCharNA <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x)
}


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
#'     Special thanks to Al J Abadi for preparing the published data in time
#'     for the 2020 BIRS Workshop, see the link here:
#'     url{https://github.com/BIRSBiointegration/Hackathon/tree/master/scNMT-seq}
#'
#' @param DataType character(1) Indicates study that produces this type of
#'     data (default: 'mouse_gastrulation')
#'
#' @param modes character() The assay types or modes of data to obtain these
#'     include single cell Chromatin Accessibilty ("acc"), Methylation ("met"),
#'     RNA-seq ("rna") by default.
#'
#' @param version character(1) The data version available in ExperimentHub
#'     defaults to the newest version ('2.0.0')
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
#' @return A single cell multi-modal \linkS4class{MultiAssayExperiment} or
#'     informative data.frame when `dry.run` is `TRUE`
#'
#' @source \url{http://ftp.ebi.ac.uk/pub/databases/scnmt_gastrulation/}
#'
#' @references
#'     Argelaguet et al. (2019)
#'
#' @examples
#' scNMT(DataType = "mouse_gastrulation", modes = "*", dry.run = TRUE)
#'
#' @export scNMT
scNMT <-
    function(
        DataType = "mouse_gastrulation", modes = "*", version = "2.0.0",
        dry.run = TRUE, verbose = TRUE, ...
    )
{
    modes_file <- system.file("extdata", "metadata.csv",
        package = "SingleCellMultiModal", mustWork = TRUE)

    DataType <- tolower(DataType)
    stopifnot(
        .isSingleCharNA(DataType), .isSingleCharNA(version)
    )

    modes_metadat <- read.csv(modes_file, stringsAsFactors = FALSE)
    filt <- modes_metadat[["DataType"]] == DataType &
        modes_metadat[["SourceVersion"]] == version
    modes_metadat <- modes_metadat[filt, ]
    eh_assays <- modes_metadat[["ResourceName"]]
    modesAvail <- .modesAvailable(eh_assays)
    resultModes <- .searchFromInputs(modes, modesAvail)
    fileIdx <- .conditionToIndex(
        resultModes, eh_assays, function(x) grepl(x, eh_assays)
    )
    fileMatches <- modes_metadat[fileIdx, c("Title", "DispatchClass")]
    eh <- .test_eh(...)

    if (dry.run) {
        return(.getResourceInfo(
            eh, modes_metadat[fileIdx, c("Title", "RDataPath")], FALSE
        ))
    }
    modes_list <- .getResources(
        eh, modes_metadat[fileIdx, c("Title", "RDataPath")], verbose
    )
    names(modes_list) <- gsub("scnmt_", "", names(modes_list))

    eh_experiments <- ExperimentList(modes_list)[resultModes]

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
