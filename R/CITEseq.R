.cord_blood <- function(ess_list)
{
    names(ess_list$experiments) <- gsub("_Counts", "", names(ess_list$experiments))
    mse <- MultiAssayExperiment::MultiAssayExperiment(experiments=(ess_list$experiments))
    return(mse)
}


#' @importFrom Matrix Matrix
.peripheral_blood <- function(ess_list)
{
    .combMatrixForAssay <- function(explist, dimslist, 
                                    assayId=c("scADT", "scHTO", "scRNA"))
    {
        match.arg(assayId)
        stopifnot( (length(grep(assayId, names(explist)))!=0), 
                   (length(grep(assayId, names(dimslist)))!=0) )
        assIdx <- grep(assayId, names(explist))
        switch(assayId,
            "scADT"=, "scHTO"={
                if(length(explist[assIdx]) == 2)
                {
                    m1 <- Matrix::Matrix(unlist(explist[assIdx]), 
                                nrow=dimslist[assIdx][[1]][1],
                                ncol=(dimslist[assIdx][[1]][2]+dimslist[assIdx][[2]][2]), 
                                sparse=TRUE)
                } else {
                    m1 <- Matrix::Matrix(explist[[assIdx]])
                }
            },
            "scRNA"={
                if(length(explist[assIdx]) == 2)
                {
                    ## we can have at last 2 matrices
                    m1 <- cbind(explist[[assIdx[1]]], explist[[assIdx[2]]])
                    
                } else {
                    m1 <- explist[[assIdx]]
                }
            },
            {stop("Unrecognized assayId: ", assayId)}
            )
        if(length(explist[assIdx]) == 2)
        {
            colnames(m1) <- c(paste0(rep(gsub("scADT|scHTO|scRNA","", 
                                            names(explist)[assIdx[1]]),
                                        dimslist[assIdx][[1]][2]), 
                                    colnames(explist[[assIdx[1]]])),
                              paste0(rep(gsub("scADT|scHTO|scRNA","", 
                                            names(explist)[assIdx[2]]),
                                        dimslist[assIdx][[2]][2]), 
                                    colnames(explist[[assIdx[2]]])))
        } else {
            colnames(m1) <- paste0(rep(gsub("scADT|scHTO|scRNA","", 
                                            names(explist)[assIdx[1]]),
                                        dimslist[assIdx][[1]][2]), 
                                    colnames(explist[[assIdx[1]]]))
        }
        return(m1)
    }
    
    .buildMap <- function(mat1, assayId)
    {
        
        map <- DataFrame(assay=assayId,
                        primary=gsub("_\\w+", "", colnames(mat1)), 
                        colname=colnames(mat1), 
                        condition=gsub("_\\w+", "", colnames(mat1)))
        return(map)
    }
    
    .buildColDat <- function(ll)
    {
        if(all(names(ll)[grep("CTCL|CTRL", names(ll))] %in% names(ll)))
        {
            cd <- data.frame(row.names=c("CTCL", "CTRL"), 
                       condition=c("Cutaneous T-cell Limphoma", "Control"))
        } else if(!isEmpty(grep("CTCL", names(ll)))) {
            cd <- data.frame(row.names=c("CTCL"), 
                              condition=c("Cutaneous T-cell Limphoma"))
        } else if(!isEmpty(grep("CTRL", names(ll)))) {
            cd <- data.frame(row.names=c("CTRL"), 
                              condition=c("Control"))
        }
        return(cd)
    }
    
    ll <- ess_list$experiments
    ll <- lapply(ll, function(x) 
    {
        x <- x[order(rownames(x)),]
    })
    dims <- lapply(ll, dim)
    expslist <- c()
    sampmap <- DataFrame()
    if(!isEmpty(grep("scADT", names(ll))))
    {
        ADTs <- .combMatrixForAssay(explist=ll, dimslist=dims, assayId="scADT")
        expslist <- c(expslist, scADT=ADTs)
        ADTsMap <- .buildMap(ADTs, assayId="scADT")
        sampmap <- rbind(sampmap, ADTsMap)
    }
    if(!isEmpty(grep("scHTO", names(ll))))
    {
        HTOs <- .combMatrixForAssay(explist=ll, dimslist=dims, assayId="scHTO")
        expslist <- c(expslist, scHTO=HTOs)
        HTOsMap <- .buildMap(HTOs, assayId="scHTO")
        sampmap <- rbind(sampmap, HTOsMap)
    }
    if(!isEmpty(grep("scRNA", names(ll))))
    {
        RNAs <- .combMatrixForAssay(explist=ll, dimslist=dims, assayId="scRNA")
        expslist <- c(expslist, scRNA=RNAs)
        RNAsMap <- .buildMap(RNAs, assayId="scRNA")
        sampmap <- rbind(sampmap, RNAsMap)
    }
    
    coldat <- .buildColDat(ll)
    mse <- MultiAssayExperiment::MultiAssayExperiment(experiments=expslist, 
                                                    sampleMap=sampmap, 
                                                    colData=coldat)
    if(!isEmpty(grep("TCR", names(ll))))
    {
        metadata(mse) <- ll[grep("TCR", names(ll))]
    }
    return(mse)
}

#' CITEseq
#' @description function assembles data on-the-fly from `ExperimentHub`
#'     to provide a \linkS4class{MultiAssayExperiment} container. Actually
#'     the `dataType` argument provides access to the available datasets
#'     associated to the package.
#' @author Dario Righelli
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
#'      \itemize{
#'         \item{peripheral_blood: } a dataset of single cells of peripheral 
#'         blood as provided in Mimitou et al. (2019).
#'         We provide two different conditions controls (CTRL) and 
#'         Cutaneous T-cell Limphoma (CTCL).
#'         Just build appropriate \code{modes} regex for subselecting the 
#'         dataset modes. 
#'          \itemize{
#'             \item{scRNA} - Mimitou scRNA-seq gene count matrix
#'             \item{scADT} - Mimitou antibody-derived tags (ADT) data
#'             \item{scHTO} - Mimitou Hashtag Oligo (HTO) data
#'             \item{TCRab} - Mimitou T-cell Receptors (TCR) alpha and beta  
#'             available through the object metadata.
#'             \item{TCRgd} - Mimitou T-cell Receptors (TCR) gamma and delta  
#'             available through the object metadata.
#'             }
#'      }
#'
#' @param DataType character(1) indicating the identifier of the dataset to
#'     retrieve.  (default "cord_blood")
#'
#' @param modes character( ) The assay types or modes of data to obtain these
#'     include scADT and scRNA-seq data by default.
#'     
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
#' @references Stoeckius et al. (2017), Mimitou et al. (2019)
#' @export
#'
#' @examples
#'
#' mse <- CITEseq(DataType="coord_blood", dry.run=FALSE)
#' experiments(mse)
#'
CITEseq <- function(DataType=c("cord_blood", "peripheral_blood"), modes="*",
                    version="1.0.0", dry.run=TRUE, verbose=TRUE, ...)
{

    ess_list <- .getResourcesList(prefix = "citeseq_", datatype = DataType,
        modes=modes, version=version, dry.run=dry.run, verbose=verbose, ...)
    if (!dry.run) {
        mse <- switch(
            DataType,
            "cord_blood" = { .cord_blood(ess_list=ess_list) },
            "peripheral_blood" = { .peripheral_blood(ess_list=ess_list) },
                ## Add here other CITE-seq datasets based on DataType identifier
            { stop("Unrecognized CITE-seq dataset name: ", DataType) }
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
#'
#' @param mse a MultiAssayExperiment object with scRNA and/or scADT and/or 
#' scHTO named experiments.
#'
#' @return a SingleCellExperiment object as widely with scRNA data as counts
#' and scADT, scHTO data as altExps.
#' If only one modality is present, it has returned as main assay of the SCE.
#' 
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment altExps
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
    stopifnot(c(is(mse, "MultiAssayExperiment"), !(length(mse)==0)))

    if(length(mse)==3)
    {
        
        scrna <- experiments(mse)[[grep("scRNA", names(mse))]]
        scadt <- SummarizedExperiment(experiments(mse)[[grep("scADT", names(mse))]])
        schto <- SummarizedExperiment(experiments(mse)[[grep("scHTO", names(mse))]])
        
        commonsamp <- intersect(intersect(colnames(scrna), colnames(scadt)), colnames(schto))
        
        schto <- schto[,(colnames(schto) %in% commonsamp)]
        scrna <- scrna[,(colnames(scrna) %in% commonsamp)]
        scadt <- scadt[,(colnames(scadt) %in% commonsamp)]
        
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts=scrna),
                                                        altExps=list(scADT=scadt, 
                                                                     scHTO=schto))
    } else if(length(mse)==2) {
        scrna <- experiments(mse)[[grep("scRNA", names(mse))]]
        if(length(grep("scADT", names(mse)))!=0)
        {
            scalt <- SummarizedExperiment(experiments(mse)[[grep("scADT", names(mse))]])
            name <- "scADT"
        } else {
            scalt <- SummarizedExperiment(experiments(mse)[[grep("scHTO", names(mse))]])
            name <- "scHTO"
        }
        commonsamp <- intersect(colnames(scrna), colnames(scalt))
        
        scalt <- scalt[,(colnames(scalt) %in% commonsamp)]
        scrna <- scrna[,(colnames(scrna) %in% commonsamp)]
        l <- list(scalt)
        names(l) <- name
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts=scrna),
                                                    altExps=l)
    } else { ## case length 1
        if(length(grep("scADT", names(mse)))!=0)
        {
            scadt <- SummarizedExperiment(experiments(mse)[[grep("scADT", names(mse))]])
            sce <- SingleCellExperiment::SingleCellExperiment(list(adt=scadt))
        } else if(length(grep("scHTO", names(mse)))!=0) {
            schto <- SummarizedExperiment(experiments(mse)[[grep("scHTO", names(mse))]])
            sce <- SingleCellExperiment::SingleCellExperiment(list(hto=schto))
        } else if(length(grep("scRNA", names(mse)))!=0) {
            scrna <- experiments(mse)[[grep("scRNA", names(mse))]]
            sce <- SingleCellExperiment::SingleCellExperiment(list(counts=scrna))
        } else if(length(grep("scHTO", names(mse)))!=0) {
            schto <- experiments(mse)[[grep("scHTO", names(mse))]]
            sce <- SingleCellExperiment::SingleCellExperiment(list(counts=schto))
        }
    }

    return(sce)
}


