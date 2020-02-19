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
uflist <- unlist(flist, use.names = FALSE)
anyDuplicated(uflist)

mets <- grepl("met", uflist)
geoaccs <- .splitselect(uflist[mets])
gsub("met|acc", "*", uflist)
filedata <- do.call(cbind.data.frame, list(geo_accession = geoaccs, split(uflist, mets),
    stringsAsFactors = FALSE))
names(filedata)[2:3] <- c("accessibilty", "methylation")
# merge(intpheno, filedata, by = "geo_accession")
## sample map for raw bsseq data
smaps <- Map(function(x, y) .nametodframe(fnames = x, gseacc = y),
    x = filedata$accessibilty, y = names(flist))
sampmap <- dplyr::bind_rows(smaps)
head(filedata)
head(sampmap)
sampmap <- sampmap[, -which(names(sampmap) %in% c("assay", "filename"))]
fmap <- merge(filedata, sampmap, by = "geo_accession")
fmap$archive <- allsups[match(fmap$GSEseries, names(allsups))]

## extract all files in tarballs
lapply(allsups[tars], function(tr) {
    untar(tr, exdir = dirname(tr))
})

access <- apply(fmap, 1L, function(x) {
    datapath <- file.path(dirname(x['archive']), x['accessibility'])
    dat <- readr::read_tsv(datapath)
    dat <- dat[complete.cases(dat), ]
    res <- GenomicRanges::makeGRangesFromDataFrame(dat,
        keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
    as(res, "GPos")
})

as(as(rlist, "GRangesList"), "RaggedExperiment")

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

