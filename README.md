
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spWAVE <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![R](https://img.shields.io/badge/R-276DC3.svg?logo=r&logoColor=white)]()
[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)]()
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**Sp**atial **W**hole **A**rea **V**ector **E**stimation (spWAVE) is
design to estimate the ligands-receptors spatial interference based on
conservative vector field model with superposition principle.

## Installation

### Dependencices

We recomannd [Seurat](https://satijalab.org/seurat/) \>= 5.0. This
package use the new Seurat ‘layers’ structure to access data, which is
not compatible with Seurat \< 5.0. For lower version Seurat and other ST
objects, we provided the function to constructe spWAVE object by
mannully giving the coordinate info, expression matrix, and other
information.

This package contains codes for Rcpp. Please ensure that the
corresponding compilation environment is installed, such as Rtools for
Windows and Xcode for macOS.

For downstream analysis, install
[tradeSeq](https://github.com/statOmics/tradeSeq) to analyze the field
quantities related differentially expressed genes.

### Install

Install from GitHub:

``` r
devtools::install_github("liuhuoz/spWAVE")
```

Or download the GitHub source zip file, and install it from local:

``` r
devtools::install_local("spWAVE.zip")
```

## Quick Start

Starting and constructing a spWAVE object from a Seurat object, you can
follow the quick start example:

``` r
library(spWAVE)
seurat_obj <- readRDS("seurat_obj.rds")
human_db <- load_database(db_source = "CellChat",db_species = "human")
#Or load a database form CellPhoneDB of mouse
#mouse_db <- load_database(db_source = "CellPhoneDB",db_species = "mouse")

spWAVE_obj <- create_spWAVE_object(seurat_obj,human_db)
spWAVE_obj %<>% perform_LR_field_calc()
spWAVE_obj %<>% perform_S2S_score_calc()
spWAVE_obj %<>% perform_C2C_score_calc(shuffle_iter=200,verbose=T)
```

## Full Tutorial for Analysis and Visualization

Please check the full tutorial in tutorial folder.

Field Estimation: “spWAVE/tutorial/Field_Est.html”

Visualization: “spWAVE/tutorial/Visualization.html”

Downstream Analysis: “spWAVE/tutorial/Downstream.html” (not write now,
will update)

All links is missing because it is not publish on Github, will be
updated in the future.
