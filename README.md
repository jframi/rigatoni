
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rigatoni

<!-- badges: start -->
<!-- badges: end -->

Rigatoni is an experimental port, based on the EBImage package, of the
grainito (<https://github.com/jframi/grainito>) ImageJ macro.

## Installation

Install dependencies from CRAN

``` r
install.packages("data.table")
install.packages("polylabelr")
install.packages("polyclip")
install.packages("devtools")
```

Install dependencies from github

``` r
devtools::install_github("MomX/Momocs2")
```

Install dependencies from Bioconductor

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EBImage")
```

Then install Rigatoni from github

``` r
devtools::install_github("jframi/rigatoni")
```

## Getting started

Read the
[Vignette](https://jframi.github.io/rigatoni/articles/Rigatoni.html) to
get started
