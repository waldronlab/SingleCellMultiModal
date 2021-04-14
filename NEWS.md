## Changes in version 1.4.0

### New features

* `SingleCellMultiModal` function allows the combination of multiple
multi-modal technologies.
* `GTseq` data from Macaulay et al. (2015) now available (@lgeistlinger)
* `SCoPE2` data from Specht et al. now available thanks to @cvanderaa (#26)
* `scMultiome` provides PBMC from 10X Genomics thanks to @rargelaguet

### Bug fixes and minor improvements

* Metadata information (function call and call to technology map) included in
`SingleCellMultiModal`
* `scNMT` includes the original call in the `MultiAssayExperiment` metadata
* Improved and edited Contributing Guidelines for clarity

## Changes in version 1.2.0

### New features

* `CITEseq` function, vignette, and 'cord_blood' data available
(@drighelli, #18)
* Include `seqFISH` function, vignette, and 'mouse_visual_cortex' data
(v1 and v2 from @drighelli, #14)
* New 'mouse_gastrulation' dataset released (version "2.0.0").
* Use `version` argument to indicate the `mouse_gastrulation` data version
* The data includes **all** cells not only the ones that passed the QC
of all three 'omics (thanks @rargelaguet, @ajabadi).

### Bug fixes and minor improvements

* Caching mechanism uses `tools::R_user_dir` and not `rappdirs`.
* Improved display of available data using `ExperimentHub` metadata.
* Improved documentation explaining versioning differences.
* Contribution guidelines available at
https://github.com/waldronlab/SingleCellMultiModal/wiki/Contributing-Guidelines
* Default `version` argument in `scNMT` function now set to "2.0.0" (version
"1.0.0" still available)

## Changes in version 1.0.0

### New features

* `scNMT` serves the mouse gastrulation dataset from Argelaguet et al. 2019
* Data set is provided by Argelaguet and colleagues via CloudStor link:
https://cloudstor.aarnet.edu.au/plus/s/Xzf5vCgAEUVgbfQ
* GitHub repository for the dataset by the authors available at:
https://github.com/rargelaguet/scnmt_gastrulation

### Bug fixes and minor improvements

* Row names in the scNMT dataset properly show mouse ENSEMBL identifiers
