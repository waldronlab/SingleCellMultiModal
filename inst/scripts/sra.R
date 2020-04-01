library(SRAdb)
library(DBI)

getGEO("GSM3440983", GSEMatrix = FALSE)
# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR8091082


library(DBI)
library(RSQLite)
library(dbplyr)
library(dplyr)

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

