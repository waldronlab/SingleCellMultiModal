# scNMT
library(GEOquery)
gse <- getGEO('GSE121708', GSEMatrix = TRUE)
# gse <- getGEO('GSE121708', GSEMatrix = FALSE)

# head(
#     pData(phenoData(gse[[1]]))
# )

# download raw data
dfdata <- getGEOSuppFiles("GSE121708")
files <- untar("GSE121708/GSE121708_RAW.tar", list = TRUE)
test <- head(files)
untar("GSE121708/GSE121708_RAW.tar", files = test)
gsm0 <- read.table("GSM3443369_E4.5-5.5_new_Plate1_A02_met.tsv.gz", header = TRUE)
gsm1 <- readr::read_tsv("GSM3443369_E4.5-5.5_new_Plate1_A02_met.tsv.gz")
gsm1.5 <- readr::read_tsv("GSM3443369_E4.5-5.5_new_Plate1_A02_acc.tsv.gz")

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
slups <- lapply(sups, getGEO, GSEMatrix = FALSE)
lapply(slups, function(x) lapply(GSMList(x), function(y) Meta(y)$title))

library(SRAdb)
# sra_con <- dbConnect(SQLite(), getSRAdbFile())
sra_con <- dbConnect(SQLite(), "SRAmetadb.sqlite")
listSRAfile("SRX4917371", sra_con, fileType = "sra")
getSRAinfo("SRX4917371", sra_con, sraType = "sra" )
getSRAfile("SRX4917371", sra_con, fileType = "sra" )
listSRAfile ("SRX4917371", sra_con, fileType = 'sra', srcType='fasp')

getSRAinfo("SRS3964271", sra_con, sraType = "sra")
getSRAinfo("SRR8091082", sra_con)
sraConvert("SRR8091082", out_type = c("study", "sample", "experiment", "run"), sra_con)

getGEO("GSM3440983", GSEMatrix = FALSE)
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8091082

metlist <- lapply(gg, getGEO)
metas <- lapply(metlist, Meta)
names(metas) <- gg
lapply(metas, function(x) as(x, "CharacterList"))


## GPL
# plat <- getGEO(annotation(gse[[1]]))


getGSEDataTables('GSE121708')

library(dbplyr)
library(dplyr)

src_dbi(sra_con)

tbl(sra_con, "run")

res <- httr::POST("https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=SRR8091082,SRR8091083")
reslist <- httr::content(res)

