# scNMT
library(GEOquery)
gse <- getGEO('GSE121708')

# head(
#     pData(phenoData(gse[[1]]))
# )

# download raw data
# dfdata <- getGEOSuppFiles("GSE121708")

# identical(
#     pData(gse[[1]])$geo_accession,
#     rownames(pData(gse[[1]]))
# )
# TRUE

gg <- head(pData(gse[[1]])$geo_accession)
metlist <- lapply(gg, getGEO)
metas <- lapply(metlist, Meta)
names(metas) <- gg
lapply(metas, function(x) as(x, "CharacterList"))


## GPL
# plat <- getGEO(annotation(gse[[1]]))


getGSEDataTables('GSE121708')
