#' Obtain a map of cell types for each dataset
#'
#' The `ontomap` function provides a mapping of all the cell names across the
#' all the data sets or for a specified data set.
#'
#' @param dataset `character()` One of the existing functions within the
#'   package. If missing, a map of all cell types in each function will
#'   be provided.
#'
#' @return A `data.frame` of metadata regarding the cell types and ontologies,
#'   if any.
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
