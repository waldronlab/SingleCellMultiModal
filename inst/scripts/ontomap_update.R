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

## reading ontology terms from kelly ontomap based on an ontomap old version
cellontokelly <- as.data.frame(
    readr::read_tsv("~/Downloads/Cell type ontology - Sheet2.tsv")
)
onto <- as.data.frame(
    readr::read_tsv("inst/extdata/ontomap.tsv")
)

## removing repetitive rows for seqFISH/scRNAseq celltypes
## aligning with newer version of ontomap
cellontokelly <-
    cellontokelly[!cellontokelly$dataset_name=="mouse_visual_cortex_scRNAseq",]
ontokey <- paste0(onto$DataType,"_",onto$function_name)
ontokey <-
    gsub("SCoPE2", "protein_SCoPE2", gsub("scMultiome", "multiome", ontokey))
ontokey <- paste0(ontokey, "_", onto$original_cell_name)
kellykey <-
    paste0(cellontokelly$dataset_name, "_", cellontokelly$original_cell_name)
## reordering
cellontokelly <- cellontokelly[match(ontokey, kellykey),]

onto$ontology_ID <- cellontokelly$ontology_ID
onto$ontology_cell_name <- cellontokelly$ontology_cell_name

## writing
write.table(
    x = onto, file = "inst/extdata/ontomap.tsv",
    quote = FALSE, sep = "\t", row.names = FALSE
)
