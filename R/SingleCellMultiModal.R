.internalMap <- data.frame(
    prefix = c("scnmt_", "pbmc_"),
    datatype = c("mouse_gastrulation", "pbmc_10x"),
    version = c("2.0.0", "1.0.0")
)

.filterMap <- function(DataTypes) {
    .internalMap[.internalMap$datatype %in% DataTypes, ]
}

#' Combining Modalities into one MultiAssayExperiment
#'
#' Combine multiple single cell modalities into one using the input of the
#' individual functions.
#'
#' @inheritParams scNMT
#'
#' @param DataTypes character() A vector of data types as indicated in each
#'     individual function by the `DataType` parameter.
#'
#' @param versions character() A vector of versions for each DataType.
#'
#' @param modes list() A list or CharacterList of modes for each data type
#'     where each element corresponds to one data type.
#'
#' @md
#'
#' @examples
#'
#' SingleCellMultiModal(c("mouse_gastrulation", "pbmc_10x"), modes = "*",
#'     version = c("2.0.0", "1.0.0"), dry.run = TRUE, verbose = TRUE
#' )
#' SingleCellMultiModal(c("mouse_gastrulation", "pbmc_10x"), modes = "*",
#'     version = c("2.0.0", "1.0.0"), dry.run = FALSE, verbose = TRUE
#' )
#'
#' @export
SingleCellMultiModal <- function(
        DataTypes, modes = "*", versions, dry.run = TRUE, verbose = TRUE, ...
    )
{
    stopifnot(is.character(DataTypes), is.character(versions))
    if (is.character(modes) && length(modes) == 1L && identical(modes, "*"))
        modes <- as.list(rep(modes, length(DataTypes)))
    resmap <- .filterMap(DataTypes)
    ess_lists <- Map(.getResourcesList, prefix = resmap[["prefix"]],
        datatype = resmap[["datatype"]], modes = modes,
        version = resmap[["version"]], dry.run = dry.run, verbose = verbose
    )
    names(ess_lists) <- resmap[["datatype"]]

    if (dry.run) { return(ess_lists) }

    new_prefix <- paste0(resmap[["datatype"]], "_")
    exps <- Map(function(x, y) {
        exs <- x[["experiments"]]
        names(exs) <- paste0(y, names(exs))
        exs
    }, x = ess_lists, y = new_prefix)
    exps <- do.call(c, unname(exps))

    samps <- Map(function(x, y) {
        samps <- x[["sampleMap"]]
        samps[["assay"]] <- paste0(y, samps[["assay"]])
        as(samps, "DataFrame")
    }, x = ess_lists, y = new_prefix)
    samps <- do.call(BiocGenerics::rbind, unname(samps))

    coldata <- Reduce(function(x, y) {
        merge(x, y, by = intersect(names(x), names(y)), all = TRUE)
    }, lapply(ess_lists, `[[`, "colData"))

    MultiAssayExperiment(
        experiments = exps,
        colData = coldata,
        sampleMap = samps
    )
}
