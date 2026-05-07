MultiEFM
========

High-dimensional multi-study elliptical factor model
=========================================================================

To robustly extract study-shared and study-specific latent features from heavy-tailed data derived from multiple heterogeneous sources, we introduce a multi-study elliptical factor model, called MultiEFM, which learns latent structures and accounts for the heterogeneity among sources without requiring finite second moments.

Installation
------------

"MultiEFM" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.

### Install from CRAN
```R
install.packages("MultiEFM")
```

### Install from github
```R
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("RRuihan-Zhang/MultiEFM")
```
Usage
-----

For usage examples and guided walkthroughs, check the `vignettes` directory of the repo or our [official website](https://rruihan-zhang.github.io/MultiEFM/).

* [Low Dimensional Example](https://rruihan-zhang.github.io/MultiEFM/articles/simu_low_dim.html)
* [High Dimensional Example](https://rruihan-zhang.github.io/MultiEFM/articles/simu_high_dim.html)

Simulated codes
---------------

For the codes in simulation study, check the `simuCodes` directory of the repo.

News
----

MultiEFM version 1.0 released! (2026-04-27)

## Authors

* **Ruihan Zhang** (Maintainer)
* **Wei Liu**
