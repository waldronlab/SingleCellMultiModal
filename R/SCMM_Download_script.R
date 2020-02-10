# scNMT
library(GEOquery)
gse <- getGEO('GSE121708', GSEMatrix = TRUE)

 head(
     pData(phenoData(gse[[1]]))
 )


allseries <- .getAllSeries(gse[[1]])

allsups <- .downloadSeries(allseries, directory = "~/data")

## checking data vs url size
vapply(allseries, .checkSize, logical(1L), directory = "~/data")

## file sizes
fs <- round(vapply(allsups, file.size, double(1L)) / 1024^3, 2)

tsvs <- grepl("tsv\\.gz$", allsups)
tars <- grepl("\\.tar$", allsups)

tt <- readr::read_tsv(allsups[6])
tp <- readr::read_tsv(allsups[7])
tr <- readr::read_tsv(allsups[1])

flist <- lapply(allsups[tars], function(x) untar(x, list = TRUE))

smaps <- Map(function(x, y) .nametodframe(fnames = x, gseacc = y),
    x = flist, y = names(flist))
sampmap <- dplyr::bind_rows(smaps)

test <- head(sampmap[sampmap$GSEseries %in% names(allsups[tars]), ])
tiny <- split(test, test$assay)

nlist <- lapply(tiny, function(fdat) {
    sampf <- setNames(fdat$filename, fdat$sample)
    lapply(sampf, readr::read_tsv)
})

onelist <- lapply(nlist, function(assay) {
    rlist <- lapply(assay, function(df) {
        dat <- head(df)
        dat <- dat[complete.cases(dat[, 1:2]), ]
        res <- GenomicRanges::makeGRangesFromDataFrame(dat,
            keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
        as(res, "GPos")
    })
    as(as(rlist, "GRangesList"), "RaggedExperiment")
})

untar("GSE121708/GSE121708_RAW.tar", files = test)

geos <- pData(gse[[1]])$geo_accession
sum(geos %in% sampmap$sample) * 100 / length(geos)

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
library(DBI)
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
library(DBI)
library(RSQLite)

sra_con <- dbConnect(SQLite(), "SRAmetadb.sqlite")
src_dbi(sra_con)

runs <- tbl(sra_con, "run")
samps <- tbl(sra_con, "sample")
group_by(runs, run_file) %>% summarize(uns = n())

sratab <- tbl(sra_con, "sra")
select(sratab, run_url_link, run_entrez_link, experiment_url_link, experiment_entrez_link, sample_url_link, sample_entrez_link, study_url_link, study_entrez_link) %>%
    collect()


res <- httr::POST("https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=SRR8091082,SRR8091083")
reslist <- httr::content(res)

