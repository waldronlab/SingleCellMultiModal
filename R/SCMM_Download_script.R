# scNMT
library(GEOquery)
gse <- getGEO('GSE121708')

# head(
#     pData(phenoData(gse[[1]]))
# )

# download raw data
dfdata <- getGEOSuppFiles("GSE121708")
files <- untar("GSE121708/GSE121708_RAW.tar", list = TRUE)

# obtain just file location
testdata <- getGEOSuppFiles("GSE121708", fetch_files = FALSE)

sampids <- vapply(strsplit(files, "_"), `[[`, character(1L), 1L)

# identical(
#     pData(gse[[1]])$geo_accession,
#     rownames(pData(gse[[1]]))
# )
# TRUE

geos <- pData(gse[[1]])$geo_accession
sum(geos %in% sampids)/length(geos) * 100

ago <- lapply(head(geos), function(acc) {
    Meta(getGEO(acc))$relation
})

notes <- notes(experimentData(gse[[1]]))$relation
sups <- gsub(".*(GSE*)", "\\1",
grep("SuperSeries", unlist(strsplit(notes, "\n")), ignore.case = TRUE, value = TRUE)
)
lapply(sups, getGEO)


library(SRAdb)
sra_con <- dbConnect(SQLite(), getSRAdbFile())
listSRAfile("SRX4917371", sra_con)

metlist <- lapply(gg, getGEO)
metas <- lapply(metlist, Meta)
names(metas) <- gg
lapply(metas, function(x) as(x, "CharacterList"))


## GPL
# plat <- getGEO(annotation(gse[[1]]))


getGSEDataTables('GSE121708')
