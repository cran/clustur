
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/SchlossLab/clustur/actions/workflows/check-full.yml/badge.svg)](https://github.com/SchlossLab/clustur/actions/workflows/check-full.yml)
[![test-coverage](https://github.com/SchlossLab/clustur/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/SchlossLab/clustur/actions/workflows/test-coverage.yml)
[![lint](https://github.com/SchlossLab/clustur/actions/workflows/lintr.yml/badge.svg)](https://github.com/SchlossLab/clustur/actions/workflows/lintr.yml)
[![pkgdown](https://github.com/SchlossLab/clustur/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/SchlossLab/clustur/actions/workflows/pkgdown.yaml)

<!-- badges: end -->

## Overview

clustur was developed for two main reasons: to allow developers to
easily create pipelines with mothur code in R, and to give users a
powerful clustering algorithm from mothur without having to download it.
Because clustur is built using rcpp, it benefits from the speed of C++
and the usability of R. Additionally, since clustur will be open source,
developers can modify are package as needed. This will not only help us
continue to create high quality packages, but will allow everyone to
take advantage of the foundation that is already built.

## Installation

For the official release from CRAN you can use the standard
`install.packages()` function:

``` r
# install via cran
install.packages("clustur")
library(clustur)
```

For the developmental version, you can use the `install_github()`
function from the {devtools} package

``` r
# install via github
devtools::install_github("SchlossLab/clustur")
library(clustur)
```

## Getting Started

To get started, look at the [“Getting
started”](https://www.schlosslab.org/clustur/) page.

## Contributions

For feature request please visit the issue page in GitHub and use the
feature tag. If you would like to contribute pull request are welcomed!
