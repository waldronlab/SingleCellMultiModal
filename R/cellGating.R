
#' addCTLabels
#'
#' @param cd the \code{colData} \code{DataFrame}
#' @param out list data structure returned by \code{getCellGroups}
#' @param outname character indicating the name of the out data structure
#' @param ct character indicating the celltype to assign in the \code{ctcol}
#' @param mkrcol character indicating the cd column to store the markers 
#' indicated by \code{outname} (default is markers)
#' @param ctcol character indicating the column in cd to store the cell type 
#' indicated by \code{ct} (default is celltype)
#' @param overwrite logical indicating if the cell types have to be overwritten 
#' without checking if detected barcodes were already assigned to other celltypes
#' @param verbose logical for having informative messages during the execution
#'
#' @return an updated version of the cd DataFrame 
#'
#' @export
addCTLabels <- function(cd, out, outname, ct, mkrcol="markers", ctcol="celltype",
                        overwrite=FALSE, verbose=TRUE)
{
    ## adds to input cd colData in the mkrcol the markers indicated by outname
    ## and in the ctcol the celltype indicated in ct
    ## the positions for the barcodes (rows in the cd) are taken in position 
    ## outname from the out structure given by function getCellGroups 
    stopifnot(any(c(mkrcol, ctcol) %in% colnames(cd)))
    stopifnot((outname %in% names(out)))
    
    cellbc <- out[[outname]]$bc
    idxc <- which(rownames(cd) %in% cellbc)
    if (length(idxc) !=0)
    {
        if (overwrite)
        {
            if(verbose) message("Blindly overwriting cell types assignements")
            cd[[mkrcol]][idxc] <- outname
            cd[[ctcol]][idxc] <- ct
        } else {
            ## checking if celltypes are already assigned
            idxnona <- which(!is.na(cd[[mkrcol]][idxc]))
            # don't get why ifelse doesn't work
            # idxcnona <- ifelse(length(idxnona)!=0, idxc[-idxnona], idxc)
            if ( length(idxnona)!=0 ) { 
                idxcnona <- idxc[-idxnona]
                if(verbose) message(length(idxnona), " Barcodes already assigned.\n",
                                    "Assigning only ", length(idxcnona), " Barcodes...")
            } else { idxcnona <- idxc }
            if (length(idxcnona)!=0) 
            {
                cd[[mkrcol]][idxcnona] <- outname
                cd[[ctcol]][idxcnona] <- ct
            } else {
                if(verbose) message("All selected Barcodes are already assigned\n",
                                    "Look at the overwrite argument to handle a more ",
                                    "brutal behaviour")
            }
        }
        
    } else {
        warning("No barcodes in cd detected for the selected ", outname, 
                "\nReturning cd as it is...")
    }
    
    return(cd)
}

#' @importFrom graphics abline smoothScatter
.plotGatingAdt <- function(mat, adt1="CD19", adt2="CD3", th1=0.2, th2=0)
{
    plot(x=mat[adt1,], y=mat[adt2,], xlab=adt1, ylab=adt2, 
         main=paste0("Gain plot with x-th: ", th1, " y-th: ", th2))
    abline(v=th1, col="red", lty=2)
    abline(h=th2, col="red", lty=2)
    smoothScatter(x=mat[adt1,], y=mat[adt2,], xlab=adt1, ylab=adt2, 
                  main=paste0("Gain plot with x-th: ", th1, " y-th: ", th2))
    
    abline(v=th1, col="red", lty=2)
    abline(h=th2, col="red", lty=2)
}


#' getCellGroups
#' 
#' @description
#' Shows the cells/barcodes in two different plots (scatter and density) 
#' divinding the space in four quadrant indicated by the two thresholds given
#' as input parameters. 
#' The x/y-axis represent respectively the two ADTs given as input.
#' It returns a list of one element for each quadrant, each with barcodes and
#' percentage (see Value section for details).
#'
#' @param mat matrix of counts or clr transformed counts for ADT data in CITEseq
#' @param adt1 character indicating the name of the marker to plot on the x-axis
#' (default is CD19).
#' @param adt2 character indicating the name of the marker to plot on the y-axis
#' (default is CD3).
#' @param th1 numeric indicating the threshold for the marker on the x-axis
#' (default is 0.2).
#' @param th2 numeric indicating the threshold for the marker on the y-axis
#' (default is 0).
#'
#' @return a list of four different element, each one indicating the quarter 
#' where the thresholds divide the plotting space, in eucledian order I, II, 
#' III, IV quadrant, indicating respectively +/+, +/-, -/+, -/- combinations 
#' for the couples of selected ADTs.
#' Each element of the list contains two objects, one with the list of detected 
#' barcodes and one indicating the percentage of barcodes falling into that
#' quadrant.
#' .
#' @details helps to do manual gating for cell type indentification with CITEseq
#' or similar data, providing cell markers. 
#' Once identified two interesting markers for a cell type, the user has to 
#' play with the thresholds to identify the cell populations specified by an 
#' uptake (+) o downtake (-) of the couple of markers (ADTs) previously selected.
#' 
#' @importFrom graphics text
#'
#' @export
getCellGroups <- function(mat, adt1="CD19", adt2="CD3", th1=0.2, th2=0)
{
    stopifnot(any(adt1,adt2) %in% rownames(mat))
    
    plot <- match.arg(plot)
    .plotGatingAdt(mat, adt1, adt2, th1, th2)
    matadt <- mat[c(adt1,adt2),]
    adt1p <- (matadt[adt1,]>th1)
    adt1m <- (matadt[adt1,]<=th1)
    adt2p <- (matadt[adt2,]>th2)
    adt2m <- (matadt[adt2,]<=th2)
    
    
    if (sum(adt1p)+sum(adt1m) != dim(mat)[2]) stop("something went wrong with adt1")
    if (sum(adt2p)+sum(adt2m) != dim(mat)[2]) stop("something went wrong with adt2")
    
    adt12pp <- which(adt1p & adt2p)
    adt12pm <- which(adt1p & adt2m)
    adt12mp <- which(adt1m & adt2p)
    adt12mm <- which(adt1m & adt2m)
    
    l <- list(
        ADT12pp=list(
            bc=colnames(matadt)[adt12pp],
            prc=((length(adt12pp)/dim(matadt)[2])*100)),
        ADT12pm=list(
            bc=colnames(matadt)[adt12pm],
            prc=((length(adt12pm)/dim(matadt)[2])*100)),
        ADT12mp=list(
            bc=colnames(matadt)[adt12mp],
            prc=((length(adt12mp)/dim(matadt)[2])*100)),
        ADT12mm=list(
            bc=colnames(matadt)[adt12mm],
            prc=((length(adt12mm)/dim(matadt)[2])*100))
    )
    names(l) <- c(paste0(adt1,"+/",adt2,"+"),
                  paste0(adt1,"+/",adt2,"-"),
                  paste0(adt1,"-/",adt2,"+"),
                  paste0(adt1,"-/",adt2,"-"))
    
    
    text((min(matadt[adt1,])+0.03), (max(matadt[adt2,])-0.05), paste(round(l[[3]]$prc), "%"))
    text((max(matadt[adt1,])-0.03), (max(matadt[adt2,])-0.05), paste(round(l[[1]]$prc), "%"))
    text((max(matadt[adt1,])-0.03), (min(matadt[adt2,])+0.05), paste(round(l[[2]]$prc), "%"))
    text((min(matadt[adt1,])+0.03), (min(matadt[adt2,])+0.05), paste(round(l[[4]]$prc), "%"))
    return(l)
}