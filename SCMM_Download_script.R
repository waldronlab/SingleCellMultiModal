#	CITE-seq: Large scale simultaneous measuremnt of epitopes and transcriptomes in single cells
library(GEOquery)
browseVignettes("GEOquery")
?getGEO
# gsm <- getGeo('GSM2695379') 
gsm <- getGEO('GSM2695379')
class(gsm)
gsm <- getGEO('GSM2695380')
