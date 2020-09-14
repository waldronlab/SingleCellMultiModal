# version 1
scmeta <- data.frame(
    DataProvider = "Dept. of Bioinformatics, The Babraham Institute, United Kingdom",
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = "https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ",
    SourceType = "RDS",
    SourceVersion = "1.0.0",
    DataType = "mouse_gastrulation",
    Maintainer  = "Marcel Ramos <marcel.ramos@roswellpark.org>",
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
    Maintainer = "Marcel Ramos <marcel.ramos@roswellpark.org>",
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
    Maintainer = c(rep("Marcel Ramos <marcel.ramos@roswellpark.org>", 2),
        rep("Dario Righelli <dario.righelli@gmail.com>", 2)),
    stringsAsFactors = FALSE
)
write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
    row.names = FALSE
)
