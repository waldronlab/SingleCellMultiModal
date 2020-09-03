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

# version 2 with spatial
scmeta <- data.frame(
    DataProvider = c(
        rep("Dept. of Bioinformatics, The Babraham Institute, United Kingdom", 2),
        "Dept. of Molecular Genetics, Allen Institute for Brain Science, United States"
    ),
    TaxonomyId = "10090",
    Species = "Mus musculus",
    SourceUrl = c(
        rep("https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ", 2),
        "https://www.dropbox.com/sh/avj4nrd4la5i88u/AACafWwBbE-xsLvOGDwRZDpYa?dl=0"
     ),
    SourceType = "RDS",
    SourceVersion = c("1.0.0", "2.0.0", "1.0.0"),
    DataType = c(rep("mouse_gastrulation", 2), "mouse_visual_cortex"),
    Maintainer = c(rep("Marcel Ramos <marcel.ramos@roswellpark.org>", 2),
        "Dario Righelli <dario.righelli@gmail.com>"),
    stringsAsFactors = FALSE
)
write.csv(
    scmeta,
    file = "inst/extdata/docuData/singlecellmultimodalv3.csv",
    row.names = FALSE
)
