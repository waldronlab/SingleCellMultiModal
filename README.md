
# SingleCellMultiModal

## Overview

`SingleCellMultiModal` is an R package that provides a convenient and
user-friendly representation of multi-modal data using
`MultiAssayExperiment`. This package introduces a suite of single-cell
multimodal landmark datasets for benchmarking and testing multimodal
analysis methods via the `ExperimentHub` Bioconductor package. The scope
of this package is to provide efficient access to a selection of
curated, pre-integrated, publicly available landmark datasets for
methods development and benchmarking.

## Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellMultiModal")
```

## Loading packages

``` r
library(SingleCellMultiModal)
library(MultiAssayExperiment)
```

# Citing SingleCellMultiModal

Your citations are crucial in keeping our software free and open source.
To cite our package see the citation (@Eckenrode2023-yq) in the
Reference section. You may also browse to the publication at [PLoS
Computational Biology](https://doi.org/10.1371/journal.pcbi.1011324).

## Representation

Users can obtain integrative representations of multiple modalities as a
`MultiAssayExperiment`, a common core Bioconductor data structure relied
on by dozens of multimodal data analysis packages.
`MultiAssayExperiment` harmonizes data management of multiple
experimental assays performed on an overlapping set of specimens.
Although originally developed for patient data from multi-omics cancer
studies, the `MultiAssayExperiment` framework naturally applies also to
single cells. A schematic of the data structure can be seen below. In
this context, “patients” are replaced by “cells”. We use
`MultiAssayExperiment` because it provides a familiar user experience by
extending `SummarizedExperiment` concepts and providing open ended
compatibility with standard data classes present in Bioconductor such as
the `SingleCellExperiment`.

<img src="https://github.com/waldronlab/MultiAssayExperiment/blob/c3c59a094e5a08111ee98b9f69579db5634d9fd4/vignettes/MultiAssayExperiment.png?raw=true" width="100%" />

# Contributions

Want to contribute to the `SingleCellMultiModal` package? We welcome
contributions from the community. Please refer to our [Contributing
Guidelines](https://github.com/waldronlab/SingleCellMultiModal/wiki/Contributing-Guidelines)
for more details.

## Further resources

For more information on the `MultiAssayExperiment` data structure,
please refer to @Ramos2017-tk as well as the [MultiAssayExperiment
vignette](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html).

# References
