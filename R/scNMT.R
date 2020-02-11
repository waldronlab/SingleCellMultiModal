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

## sample map for raw bsseq data
smaps <- Map(function(x, y) .nametodframe(fnames = x, gseacc = y),
    x = flist, y = names(flist))
sampmap <- dplyr::bind_rows(smaps)

all(unlist(flist, use.names = FALSE) %in% sampmap$filename)

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

## sample map for RNA seq data
allpheno <- lapply(allseries, function(x) {
    gse <- getGEO(x, GSEMatrix = TRUE)
    pdat <- pData(phenoData(gse[[1]]))
    pdat$ex.plate <- gsub("Sample [0-9]*_", "",
        gsub("\\s+\\(sc[A-Z]+-Seq\\)", "", pdat$title))
    pdat$geo_series <- x
    pdat
})

cnames <- Reduce(intersect, lapply(allpheno, names))
intpheno <- lapply(allpheno, `[`, cnames)

intpheno <- dplyr::bind_rows(intpheno)
