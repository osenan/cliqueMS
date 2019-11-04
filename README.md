
Introduction
============

"cliqueMS" annotates processed LC/MZ data. This R package obtains annotation for isotopes, ion adducts and fragmentation adducts. The adducts list can be supplied by the user or either use one of the package's lists.

Installation
============

Installation has been tested on Linux, Windows and macOS.

You can install "cliqueMS" 'devel' version from github with:

``` r
# install.packages("devtools")
devtools::install_github("osenan/cliqueMS")
```

Also, you can install it the release version from Bioconductor with:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("cliqueMS")
```

Usage and examples
==================

For a tutorial on "cliqueMS" please check the vignette:

``` r
browseVignettes("cliqueMS")
```
