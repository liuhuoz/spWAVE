
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spWAVE <img src="../man/figures/logo.png" align="right" height="138" alt="" />

<div align="left">

<a href="../README.md">English</a>  \|  <a href="README_CH.md">中文</a>

</div>

<!-- badges: start -->

[![R](https://img.shields.io/badge/R-276DC3.svg?logo=r&logoColor=white)]()
[![stability-beta](https://img.shields.io/badge/stability-beta-33bbff.svg)]()
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

spWAVE(**sp**atial **W**hole **A**rea **V**ector **E**stimation)
是基于保守力场模型与叠加原理所开发的一个空间转录组分析工具，
用于估计配体-受体在空间中的相互作用场。

<div align="center">

<img src="../man/figures/sample1.png" width="80%"></img>

</div>

## 安装

### 依赖

我们推荐 [Seurat](https://satijalab.org/seurat/) 版本大于 5.0。 由于
spWAVE 通过新的Seurat的 “layer” 数据结构来获取数据，因而不兼容低于 5.0
的 Seurat 对象。 对于低版本的 Seurat
和其他空间转录组存储结构，我们也开放了输入参数来帮助构筑spWAVE对象。
通过手动输入符合格式要求的空间坐标、表达矩阵和其他信息，同样可以轻松构筑spWAVE对象。

由于本软件包包含 Rcpp 的代码，因此请确保具有相应的编译环境，如 Windows
的 Rtools 和 macOS 的 Xcode。

我们提供 [tradeSeq](https://github.com/statOmics/tradeSeq)
分析接口，用于场物理量相关的差异表达基因。 如有需要，请先安装
[tradeSeq](https://github.com/statOmics/tradeSeq) 。

### 安装

从 GitHub 安装:

``` r
devtools::install_github("liuhuoz/spWAVE")
```

或是下载 GitHub 源代码压缩包进行本地安装:

``` r
devtools::install_local("spWAVE.zip")
```

## 快速上手

对于细胞数小于10000的数据集， 可以通过 Seurat V5 对象快速构建 spWAVE
对象并进行分析， 依照如下代码：

``` r
library(spWAVE)
seurat_obj <- readRDS("seurat_obj.rds")
human_db <- load_database(db_source = "CellChat",db_species = "human")
#或是选择CellPhoneDB中的数据库，如小鼠的数据库：
#mouse_db <- load_database(db_source = "CellPhoneDB",db_species = "mouse")

spWAVE_obj <- create_spWAVE_object(seurat_obj,human_db)
spWAVE_obj %<>% perform_LR_field_calc()
spWAVE_obj %<>% perform_S2S_score_calc()
spWAVE_obj %<>% perform_C2C_score_calc(shuffle_iter=200,verbose=T)
```

对于细胞数大于10000的数据集，我们推荐使用“挖孔”法进行计算，减少计算资源的消耗。

``` r
spWAVE_obj <- create_spWAVE_object(seurat_obj,human_db)

spWAVE_obj %<>% generate_kmeans_coord()
spWAVE_obj %<>% generate_meta_expr()

spWAVE_obj %<>% perform_LR_field_hole_calc()
spWAVE_obj %<>% perform_S2S_score_calc()
spWAVE_obj %<>% perform_C2C_score_calc(shuffle_iter=200,verbose=T)
```

## 分析与可视化的完整教程

请查看tutorial文件夹下的相应的完整教程

Please check the full tutorial in tutorial folder.

低分辨率空转数据集的分析: [Analysis LR field of 10x Visium Mouse Brain
dataset using spWAVE](../tutorial/Field_Est2.html)

高分辨率空转数据集分析与进阶技巧: [Analysis LR field of large dataset
using spWAVE](../tutorial/Field_Est3.html)

下游分析: “spWAVE/tutorial/Downstream.html” (待完成)
