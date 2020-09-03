## Changes in version 1.2.0

### New features

* New 'mouse_gastrulation' dataset released (version "2.0.0").
* Use `version` argument to indicate the `mouse_gastrulation` data version
* The data includes **all** cells not only the ones that passed the QC
of all three 'omics (thanks @rargelaguet, @ajabadi).

### Bug fixes and minor improvements

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
