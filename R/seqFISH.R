
#' seqFISH
#' 
#' @description seqFISH function assembles data on-the-fly from `ExperimentHub` 
#'     to provide a \linkS4class{MultiAssayExperiment} container. Actually 
#'     the `dataType` argument provides access to the available datasets 
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
#'         \itemize{
#'             \item{scRNA_Counts} - Tasic scRNA-seq gene count matrix
#'             \item{scRNA_Labels} - Tasic scRNA-seq cell labels
#'             \item{seqFISH_Coordinates} - Zhu seq-FISH spatial coordinates
#'             \item{seqFISH_Counts} - Zhu seq-FISH gene counts matrix
#'             \item{seqFISH_Labels} - Zhu seq-FISH cell labels
#'         }
#'     }
#'     
#' @param DataType character( ) indicating the identifier of the dataset to 
#'     retrieve.  (default "mouse_visual_cortex")
#'
#' @param modes character( ) The assay types or modes of data to obtain these
#'     include seq-FISH and scRNA-seq data by default.
#'
#' @param version character(1) The data version available in ExperimentHub
#'     defaults to the newest version.
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
#' @return A MultiAssayExperiment class object with seq-FISH data loaded
#' 
#' @author Dario Righelli <dario.righelli <at> gmail.com>
#' @importFrom SpatialExperiment SpatialExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
#' @export
#' @examples
#' seqFISH(dataType = "mouse_visual_cortex", modes = "*", dry.run = TRUE)
seqFISH <- 
    function(
        DataType="mouse_visual_cortex", modes="*", version, 
        dry.run=TRUE, verbose=TRUE, ...
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
    if (identical(modes, "*") && dry.run) {
        message("Available data modes for\n",
                "  ", DataType, ":\n",
                paste(
                    strwrap(paste(modesAvail, collapse = ", "),
                            width = 46, indent = 4, exdent = 4),
                    collapse = "\n"
                )
        )
        return(invisible())
    }
    
    
    resultModes <- .searchFromInputs(modes, modesAvail)
    fileIdx <- .conditionToIndex(resultModes, eh_assays, 
                                function(x) grepl(x, eh_assays))
    fileMatches <- modes_metadat[fileIdx, c("Title", "DispatchClass")]
    
    
    if (dry.run) { return(fileMatches) }
    eh <- .test_eh(...)
    modes_list <- .getResources(
        eh, modes_metadat[fileIdx, c("Title", "RDataPath")], verbose
    )
    
    names(modes_list) <- gsub("seqfish_", "", names(modes_list))
    
    ess_names <- c("colData", "metadata", "sampleMap")
    
    ess_idx <- .conditionToIndex(ess_names, eh_assays,
                                 function(x) grepl(x, eh_assays))
    
    ess_list <- .getResources(eh,
                    modes_metadat[ess_idx, c("Title", "RDataPath")], verbose)
    
    names(ess_list) <- gsub("seqfish_", "", names(ess_list))

    switch (DataType,
        "mouse_visual_cortex" = {
            mse <- .mouse_visual_cortex(modes_list=modes_list,
                                        version=version)
        },
        ## Add here other seqFISH datasets based on DataType identifier
        {
            stop("Unrecognized seqFISH dataset name!")
        }
    )

    return(mse)
}

.mouse_visual_cortex <- function(modes_list, version)
{
    switch(version,
        "1.0.0" = {
            sce <- SingleCellExperiment::SingleCellExperiment(
                rowData=rownames(modes_list$scRNA_Full_Counts),
                colData=modes_list$scRNA_Full_Labels,
                assays=S4Vectors::SimpleList(
                    counts=as.matrix(modes_list$scRNA_Full_Counts)))
        },
        "2.0.0" = {
            sce <- SingleCellExperiment::SingleCellExperiment(
                rowData=rownames(modes_list$scRNA_Counts),
                colData=modes_list$scRNA_Labels,
                assays=S4Vectors::SimpleList(
                    counts=as.matrix(modes_list$scRNA_Counts)))
        }
    )
    se <- SpatialExperiment::SpatialExperiment(
        rowData=rownames(modes_list$seqFISH_Counts),
        colData=modes_list$seqFISH_Labels,
        assays=S4Vectors::SimpleList(
            counts=as.matrix(modes_list$seqFISH_Counts)),
        spatialCoords=modes_list$seqFISH_Coordinates)
    
    
    
    mse <- MultiAssayExperiment::MultiAssayExperiment(
        experiments=c("seqFISH"=se, "scRNAseq"=sce))
    return(mse)
}