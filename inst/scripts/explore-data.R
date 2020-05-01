.checkSize <- function(access, directory = ".", verbose = TRUE) {
    finfo <- getGEOSuppFiles(access, fetch_files = FALSE)
    fur <- as.character(finfo[['url']])

    header <- httr::HEAD(fur)$headers
    header_bytes <- as.numeric(header$`content-length`)

    fpath <- file.path(directory, access, finfo[['fname']])
    local_bytes <- file.size(fpath)

    if (verbose)
        message("url: ", header_bytes, " vs. local: ", local_bytes)
    identical(header_bytes, local_bytes)
}

.getAllSeries <- function(gseacc) {
    if (!is(gseacc, "ExpressionSet"))
        stop("Provide an 'ExpressionSet' for obtaining series from 'colnames'")
    cells <- colnames(gseacc)
    metalist <- lapply(cells, function(x)
        Meta(getGEO(x, GSEMatrix = FALSE))
    )
    allseries <- unique(unlist(
        lapply(metalist, `[[`, "series_id")
    ))
}

.downloadSeries <- function(series, directory = ".") {
    directory <- normalizePath(directory)
    vapply(series, function(x) {
        finfo <- getGEOSuppFiles(x, fetch_files = FALSE)
        fname <- finfo[['fname']]
        fpath <- file.path(directory, x, fname)
        if (file.exists(fpath) && .checkSize(x, directory, verbose = FALSE))
            message("File exists: ", basename(fpath))
        else
            getGEOSuppFiles(x, baseDir = directory)
        fpath
    }, character(1L))
}

.nametodframe <- function(fnames, splitter = "_") {
    filen <- fnames
    fnames <- strsplit(fnames, splitter)
    datf <- lapply(fnames, function(fname) {
        fname[length(fname)] <-
            gsub("\\.tsv\\.gz", "", fname[length(fname)])
        samp <- fname[[1L]]
        inst <- fname %in% c("new", "repeat")
        instance <- ifelse(any(inst), fname[inst], NA_character_)
        fname <- fname[!inst]
        plates <- c(1L:(length(fname)-3L), length(fname))
        expro <- c(1L, (length(fname)-2):length(fname))
        explay <- c(1L, length(fname))

        plate <- paste0(fname[-plates], collapse = "_")
        extract <- paste0(fname[-expro], collapse = "_")
        center <- paste0(fname[-explay], collapse = "_")
        c(samp, extract, instance, plate, center)
    })
    sampmap <- do.call(rbind.data.frame, c(datf, stringsAsFactors = FALSE))
    names(sampmap) <- c("geo_accession", "extraction", "instance", "plate",
        "ext_plate")
    sampmap
}

.splitselect <- function(char, indx = 1L, sep = "_") {
    vapply(strsplit(char, "_"), `[[`, character(1L), indx)
}

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
# lapply(allsups[tars], function(tr) {
#     untar(tr, exdir = dirname(tr))
# })
.load_raw_file <- function(filepath)
{
    ctypes <- list(chr = col_character(), pos = col_double(),
        met_reads = col_double(), nonmet_reads = col_double(),
        rate = col_double())
    readr::read_tsv(filepath, col_types = ctypes)
}

projpath <- file.path("~/data", "scNMT")
dir.create(projpath, recursive = TRUE)

library(readr)
access <- apply(fmap[1:5, ], 1L, function(x) {
    datapath <- file.path("~/data", x['GSEseries'], x['accessibility'])
    stopifnot(file.exists(datapath))
    dat <- .load_raw_file(datapath)
    if (any(is.na(dat[, "pos"])))
        warning("File contains missing: ", datapath)
    dat <- dat[complete.cases(dat), ]
#    dtype <- file.path(projpath, "ACC")
#    dir.create(dtype, recursive = TRUE)
#    bychr <- split(dat, dat[["chr"]])
#    Map(function(cherf, chername) {
#        readr::write_tsv(cherf, file.path(dtype, paste0(chername, "_pos.txt")))
#    }, cherf = bychr,  chername = unique(dat[["chr"]]))
    res <- GenomicRanges::makeGRangesFromDataFrame(dat,
        keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")
    as(res, "GPos")
})

names(access) <- fmap[["ext_plate"]][1:5]
grltest <- as(access, "GRangesList")
library(RaggedExperiment)
accmet <- as(grltest, "RaggedExperiment")

system.time({
  ca <- compactAssay(accmet, "rate")
})

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

intpheno[intpheno$geo_series == "GSE133689", "ext_plate"] <-
    paste0("TET_", intpheno[intpheno$geo_series == "GSE133689", "ext_plate"])


## match GSM names with plates in rna data
seriecols <- lapply(setNames(allsubseries[tsvs], allsubseries[tsvs]),
    function(serie) {
        platematchidx <-
            match(colnames(rnalist[[serie]])[-1], intpheno[["ext_plate"]])
        intpheno[ platematchidx, "geo_accession" ]
    }
)
lapply(seriecols, head)

rnalist2 <- Map(function(x, y) {
    names(x)[-1] <- y
    x
}, x = rnalist, y = seriecols)

head(names(rnalist2[[2]]), 20)

tt <- unname(split(t(combn(1:3, 2)), rep(1:2, each = 3)))

mapply(function(x, y) {
    sum(names(rnalist[[x]])[-1] %in% names(rnalist[[y]])[-1])
}, x = tt[[2]], y = tt[[1]])

ext <- ExperimentList(
    lapply(rnalist2, function(x) {
        rnames <- x[[1]]
        dat <- data.matrix(x[, -1])
        rownames(dat) <- rnames
        SummarizedExperiment(assays = list(counts = dat))
    })
)

rownames(intpheno) <- intpheno$geo_accession

MultiAssayExperiment(ext, intpheno)

