# scNMT
library(GEOquery)
gse <- getGEO('GSE121708', GSEMatrix = TRUE)

 head(
     pData(phenoData(gse[[1]]))
 )


# Use metadata to search for all series
# allseries <- .getAllSeries(gse[[1]])
allseries <- c("GSE121650", "GSE121708", "GSE121690",
    "GSE133687", "GSE133688", "GSE133689", "GSE133725")
allsubseries <- c("GSE121650", "GSE121690", "GSE133687",
    "GSE133688", "GSE133689", "GSE133725")

source("utils.R")

# Efficient download supplemental data and get path
allsups <- .downloadSeries(allsubseries, directory = "~/data")

## checking data vs url size
vapply(allseries, .checkSize, logical(1L), directory = "~/data")

## file sizes
round(vapply(allsups, file.size, double(1L)) / 1024^3, 2)

tsvs <- grepl("tsv\\.gz$", allsups)
tars <- grepl("\\.tar$", allsups)

rnalist <- lapply(allsups[tsvs], readr::read_tsv)

flist <- lapply(allsups[tars], function(x) untar(x, list = TRUE))
lapply(flist, head)
flist <- stack(flist)
names(flist) <- c("filename", "GSEseries")
flist[] <- lapply(flist, as.character)
head(flist)

## check duplicates
anyDuplicated(flist$filename)

mets <- grepl("met", flist$filename)
flist <- do.call(cbind.data.frame, split(flist, mets))
flist <- flist[, -which(names(flist) == "FALSE.GSEseries")]
names(flist) <- c("accessibility", "methylation", "GSEseries")
flist$geo_accession <- .splitselect(flist$accessibility)
head(flist)

dfromfile <- .nametodframe(flist$accessibility)
fmap <- merge(flist, dfromfile, by = "geo_accession")

## samplemap
head(fmap)

## extract all files in tarballs
lapply(allsups[tars], function(tr) {
    untar(tr, exdir = dirname(tr))
})

access <- apply(fmap, 1L, function(x) {
    datapath <- file.path("~/data", x['GSEseries'], x['accessibility'])
    stopifnot(file.exists(datapath))
    dat <- readr::read_tsv(datapath)
    dat <- dat[complete.cases(dat), ]
    res <- GenomicRanges::makeGRangesFromDataFrame(dat,
        keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
    as(res, "GPos")
})
names(access) <- fmap[["ext_plate"]]
grltest <- as(access, "GRangesList")
accmet <- as(grltest, "RaggedExperiment")


allGEOs <- lapply(allsubseries, getGEO, GSEMatrix = TRUE)
allGEOs <- unlist(allGEOs, recursive = TRUE)
names(allGEOs) <- allsubseries
## sample map for RNA seq data
allpheno <- lapply(allsubseries, function(x) {
    gse <- getGEO(x, GSEMatrix = TRUE)
    pdat <- pData(phenoData(gse[[1]]))
    pdat$ext_plate <- gsub("Sample [0-9]*_", "",
        gsub("\\s+\\(sc[A-Z]+-Seq\\)", "", pdat$title))
    pdat$geo_series <- x
    pdat
})

cnames <- Reduce(intersect, lapply(allpheno, names))
intpheno <- lapply(allpheno, `[`, cnames)

intpheno <- dplyr::bind_rows(intpheno)

lapply(rnalist, function(x)
    all(names(x[-1]) %in% intpheno$ext_plate)
)

lapply(rnalist, function(x) anyDuplicated(names(x)))

tt <- unname(split(t(combn(1:3, 2)), rep(1:2, each = 3)))
mapply(function(x, y) {
    sum(names(rnalist[[x]])[-1] %in% names(rnalist[[y]])[-1])
}, x = tt[[2]], y = tt[[1]])

ext <- ExperimentList(
    lapply(rnalist, function(x) {
        rnames <- x[[1]]
        dat <- data.matrix(x[, -1])
        rownames(dat) <- rnames
        SummarizedExperiment(assays = list(counts = dat))
    })
)

