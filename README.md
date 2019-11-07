
Introduction
============

"cliqueMS" annotates processed LC/MZ data. This R package obtains annotation for isotopes, ion adducts and fragmentation adducts. The adducts list can be supplied by the user or either use one of the package's lists.

Installation
============

Installation has been tested on Linux, Windows and macOS.

You can install it the "release" version from Bioconductor with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cliqueMS", version = "3.10")
```
To install the "devel" version just change the `version` parameter to 3.11.

Alternatively, you can install "cliqueMS" "devel" version from github with:

``` r
# install.packages("devtools")
devtools::install_github("osenan/cliqueMS")
```


Usage and examples
==================

For a tutorial on "cliqueMS" please check the vignette:

``` r
browseVignettes("cliqueMS")
```
