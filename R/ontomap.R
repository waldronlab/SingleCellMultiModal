#' Obtain a map of cell types for each dataset
#'
#' The `ontomap` function provides a mapping of all the cell names across the
#' all the data sets or for a specified data set.
#'
#' @param dataset `character()` One of the existing functions within the
#'   package. If missing, a map of all cell types in each function will
#'   be provided.
#'
#' @details
#' Note that `CITEseq` does not have any cell annotations; therefore, no entries
#' are present in the `ontomap`.
#'
#' @return A `data.frame` of metadata with cell types and ontologies
#'
#' @examples
#'
#' ontomap(dataset = "scNMT")
#'
#' @export
ontomap <- function(
    dataset = c("scNMT", "scMultiome", "SCoPE2", "CITEseq", "seqFISH")
) {
    dataset <- match.arg(dataset, several.ok = TRUE)
    omap <- system.file(
        "extdata", "ontomap.tsv",
        package = "SingleCellMultiModal", mustWork = TRUE
    )
    map <- utils::read.delim(omap)
    dnames <- map[["function_name"]]
    map[dnames %in% dataset, ]
}
