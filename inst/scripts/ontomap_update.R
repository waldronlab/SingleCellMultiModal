## read in
onto <- readr::read_tsv("inst/extdata/ontomap.tsv")

## modification
onto <- as.data.frame(onto)
onto[onto$DataType == "macrophage_differentiation_protein", "DataType"] <-
    "macrophage_differentiation"

## output checking
stopifnot(
    identical(length(unique(onto[["DataType"]])), 4L)
)

## writing
write.table(
    x = onto, file = "inst/extdata/ontomap.tsv",
    quote = FALSE, sep = "\t", row.names = FALSE
)
