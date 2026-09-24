#This file is used to declare global variables used in spWAVE
#Mainly used in functions of ggplot2 and dplyr because of the tidyverse system.
utils::globalVariables(
  c("x","y","Ex","Ex.u","Ey","Ey.u",
    "Family","Ligand","Receptor","S2T","Source","Target",
    "lig","rec","p_value",
    "Type2","U","q_net","raw_score","scale_score",
    "waldStat",
    "barcode","position","spot_id","uni_id",
    "xend", "yend"
    )
)

#This part is used for handle the Rcpp related functions
#' @useDynLib spWAVE, .registration = TRUE
#' @exportPattern "^[[:alpha:]]+"
#' @importFrom Rcpp sourceCpp
NULL
