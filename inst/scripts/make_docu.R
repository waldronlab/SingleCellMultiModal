# version 1
scmeta <- data.frame(
    DataProvider = "Dept. of Bioinformatics, The Babraham Institute, United Kingdom",
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = "https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ",
    SourceType = "RDS",
    SourceVersion = "1.0.0",
    DataType = "mouse_gastrulation",
    Maintainer  = "Ricard Argelaguet <ricard@ebi.ac.uk>",
    stringsAsFactors = FALSE
)
write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv1.csv",
    row.names = FALSE
)

# version 2
scmeta <- data.frame(
    DataProvider =
        "Dept. of Bioinformatics, The Babraham Institute, United Kingdom",
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = "https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ",
    SourceType = "RDS",
    SourceVersion = c("1.0.0", "2.0.0"),
    DataType = "mouse_gastrulation",
    Maintainer  = "Ricard Argelaguet <ricard@ebi.ac.uk>",
    stringsAsFactors = FALSE
)
write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv2.csv",
    row.names = FALSE
)

# version 3 with spatial
scmeta <- data.frame(
    DataProvider = c(
        rep("Dept. of Bioinformatics, The Babraham Institute, United Kingdom", 2),
        rep("Dept. of Molecular Genetics, Allen Institute for Brain Science, United States", 2)
    ),
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = c(
        rep("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ", 2),
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585",
        "https://www.dropbox.com/sh/avj4nrd4la5i88u/AACafWwBbE-xsLvOGDwRZDpYa?dl=0"
    ),
    SourceType = c("RDS", "RDS", "TXT", "TXT"),
    SourceVersion = c("1.0.0", "2.0.0", "1.0.0", "2.0.0"),
    DataType = c(rep("mouse_gastrulation", 2), rep("mouse_visual_cortex", 2)),
    Maintainer = c(rep("Ricard Argelaguet <ricard@ebi.ac.uk>", 2),
                   rep("Dario Righelli <dario.righelli@gmail.com>", 2)),
    stringsAsFactors = FALSE
)
write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
    row.names = FALSE
)


# version 4 with cord_blood
scmeta <- data.frame(
    DataProvider = c(
        rep("Dept. of Bioinformatics, The Babraham Institute, United Kingdom", 2),
        rep("Dept. of Molecular Genetics, Allen Institute for Brain Science, United States", 2),
        "Innovation Lab, New York Genome Center, New York, United States"
    ),
    TaxonomyId = c(rep("10090",4), "9606"),
    Species = c(rep("Mus musculus", 4), "Homo sapiens"),
    SourceUrl = c(
        rep("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ", 2),
        rep("https://www.dropbox.com/sh/avj4nrd4la5i88u/AACafWwBbE-xsLvOGDwRZDpYa?dl=0", 2),
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866"
    ),
    SourceType = c(rep("RDS", 2), rep("TXT",3)),
    SourceVersion = c("1.0.0", "2.0.0", "1.0.0", "2.0.0", "1.0.0"),
    DataType = c(rep("mouse_gastrulation", 2), rep("mouse_visual_cortex",2), "coord_blood"),
    Maintainer = c(rep("Ricard Argelaguet <ricard@ebi.ac.uk>", 2),
                   rep("Dario Righelli <dario.righelli@gmail.com>",3)),
    stringsAsFactors = FALSE
)

write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
    row.names = FALSE
)

# indv cord_blood
citeseqmeta <- data.frame(
    DataProvider =
        "Innovation Lab, New York Genome Center, New York, United States",
    TaxonomyId = "9606",
    Species = "Homo sapiens",
    SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866",
    SourceType = "TXT",
    SourceVersion = "1.0.0",
    DataType = "coord_blood",
    Maintainer = "Dario Righelli <dario.righelli@gmail.com>",
    stringsAsFactors = FALSE
)

write.csv(
    citeseqmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv5.csv",
    row.names = FALSE
)
#
#
# # version 2 with spatial
# scmeta <- data.frame(
#     DataProvider = c(
#         rep("Dept. of Bioinformatics, The Babraham Institute, United Kingdom", 2),
#         rep("Dept. of Molecular Genetics, Allen Institute for Brain Science, United States", 2),
#         "Innovation Lab, New York Genome Center, New York, United States"
#     ),
#     TaxonomyId = c(rep("10090",4), "9606"),
#     Species = c(rep("Mus musculus", 4), "Homo sapiens"),
#     SourceUrl = c(
#         rep("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ", 2),
#         rep("https://www.dropbox.com/sh/avj4nrd4la5i88u/AACafWwBbE-xsLvOGDwRZDpYa?dl=0", 2),
#         "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866"
#     ),
#     SourceType = c(rep("RDS", 2), rep("TXT",3)),
#     SourceVersion = c("1.0.0", "2.0.0", "1.0.0", "2.0.0", "1.0.0"),
#     DataType = c(rep("mouse_gastrulation", 2), rep("mouse_visual_cortex",2), "coord_blood"),
#     Maintainer = c(rep("Marcel Ramos <marcel.ramos@roswellpark.org>", 2),
#                    rep("Dario Righelli <dario.righelli@gmail.com>",3)),
#     stringsAsFactors = FALSE
# )
# write.csv(
#     scmeta,
#     file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
#     row.names = FALSE
# )

# version 5 pbmc
scmeta <- data.frame(
    DataProvider = "European Bioinformatics Institute (EMBL-EBI), United Kingdom",
    TaxonomyId = "9606",
    Species = "Homo sapiens",
    SourceUrl = "http://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/filtered_feature_bc_matrix/",
    SourceVersion = "1.0.0",
    DataType = "pbmc_10x",
    Maintainer  = "Ricard Argelaguet <ricard@ebi.ac.uk>",
    stringsAsFactors = FALSE
)

write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv6.csv",
    row.names = FALSE
)

## version 7: creating metadata for the SCoPE2 dataset
scope2meta <- data.frame(
    DataProvider = paste0("Slavov Laboratory and SCP Center at ",
                          "Northeastern University, Boston, United ",
                          "states"),
    TaxonomyId = "9606",
    Species = "Homo sapiens",
    SourceUrl = c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142392",
                  "https://drive.google.com/file/d/1sF5STkofF_f2msnYaaYdWabou84Qf2Xr/view?usp=sharing",
                  "https://drive.google.com/file/d/16vf6rjIsk-oK9naAH6BQnCFrlWnYtJsS/view?usp=sharing"),
    SourceType = c("CSV", "CSV", "CSV"),
    SourceVersion = "1.0.0",
    DataType = "macrophage_differentiation",
    Maintainer = "Christophe Vanderaa <christophe.vanderaa@uclouvain.be>",
    stringsAsFactors = FALSE
)

write.csv(
    scope2meta,
    file = "inst/extdata/docuData/singlecellmultimodalv7.csv",
    row.names = FALSE
)
