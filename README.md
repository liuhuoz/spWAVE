
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spWAVE <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![R](https://img.shields.io/badge/R-276DC3.svg?logo=r&logoColor=white)]()
[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)]()
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Sp**atial **W**hole **A**rea **V**ector **E**stimation (spWAVE) is
design to estimate the ligand receptor spatial interference based on
vector field superposition principle.

## Installation

### Dependency

Require [Seurat](https://satijalab.org/seurat/) \>= 5.0, because the
data structure changed after 5.0, This package use the new ‘layers’
structure to access data, which is not compatible with Seurat \< 5.0.
And we recomannd Seurat \>=5.1.0, due to a better compatible with new ST
techniques, escpecially for the HD single-cell ST.

### Install

Right now, it is a non-open source package, even not published on
github, so please download the .zip file and install it from local.

``` r
devtools::install_local("spWAVE.zip")
```

## Tutorial

Please check the tutorial in this package folder.

Field Estimation: “spWAVE/tutorial/Field_Est.html”

Visualization: “spWAVE/tutorial/Visualization.html”

Downstream Analysis: “spWAVE/tutorial/Downstream.html” (not write now,
will update)

All links is missing because it is not publish on Github, will update in
the future.
