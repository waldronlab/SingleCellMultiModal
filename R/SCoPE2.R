SCoPE2 <- function(DataType="SCoPE2", 
                   modes="*", 
                   version,
                   dry.run=TRUE, 
                   verbose=TRUE, 
                   ...) {
    ## Retrieve the different resources from ExperimentHub
    ess_list <- .getResourcesList(prefix = "SCoPE2_", 
                                  datatype = DataType,
                                  modes = modes, 
                                  version = version, 
                                  dry.run = dry.run,
                                  verbose = verbose, 
                                  ...)
    ## If dry.run, return only the information table
    if (dry.run) return(ess_list)
    ## The retrieved data is stored in the `experiments` element
    modes_list <- ess_list[["experiments"]]
    ## Get a formated MultiAssayExperiment object from the list of 
    ## data resources
    .SCoPE2(modes_list = modes_list)
}

## Internal function that builds an MAE object from the data pieces 
## distributed on EH. 
.SCoPE2 <- function(modes_list = modes_list) {
    ## Construct the scRNA-Seq data
    ## Batch 1
    rna1 <- HDF5Array(file = modes_list[[1]], name = "rna1")
    sc1 <- SingleCellExperiment(assays = list(counts = rna1))
    ## Batch 2
    rna2 <- HDF5Array(file = modes_list[[1]], name = "rna2")
    sc2 <- SingleCellExperiment(assays = list(counts = rna1))
    ## Construct the SCP data
    scp <- SingleCellExperiment(assays = modes_list[[2]],
                                colData = modes_list[[3]])
    ## Build the MAE
    MultiAssayExperiment(experiments = ExperimentList(scRNAseq1 = sc1,
                                                      scRNAseq2 = sc2,
                                                      scProtein = scp))
}