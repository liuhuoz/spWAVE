
#' The spWAVE Class
#'
#' The spWAVE object is created from a spatial transcriptomic Seurat V5 object.
#' When inputting an data matrix, it takes a digital data matrices as input. Genes should be in rows and cells in columns. rownames and colnames should be included.
#' The class provides functions for data preprocessing, intercellular communication network inference, communication network analysis, and visualization.
#'
#'
#'# Class definitions
#' @importFrom methods setClassUnion
#' @importClassesFrom Matrix dgCMatrix
setClassUnion(name = 'Mat_like', members = c("matrix", "dgCMatrix","data.frame"))


#' slots of class spWAVE
#'
#' @slot expr_raw unfiltered normalized expr matrix for seurat
#' @slot expr_complex calculated and filtered receptors complex expr matrix
#' @slot kept_db database after filtering
#' @slot coord spatial coordinates of cells
#' @slot cluster_info functional annotation cluster info of cells
#' @slot meta_coord geometric center of meta clustered cells' coordinates 
#' @slot meta_coord_clu cells barcode with meta cluster info
#' @slot single_mole_field single molecular field list
#' @slot LR_pair_field ligand-receptor pair field list
#' @slot LR_family_field ligand-receptor family field list
#' 
#' @exportClass spWAVE
#' @importFrom methods setClass
spWAVE <- 
  methods::setClass("spWAVE", slots = 
    c(
      #** setup
      expr_raw = 'Mat_like',
      expr_complex = 'Mat_like',
      kept_db = "data.frame",
      coord = "data.frame",
      cluster_info = "data.frame",
      #** pre-process
      meta_coord = "data.frame",
      meta_coord_clu = "data.frame",
      #meta_complex = "Mat_like",
      #** field result
      single_mole_field = "list",
      LR_pair_field = "list",
      LR_family_field = "list"
    )
)


#' show method for spWAVE
#'
#' @param spWAVE object
#' @param show show the object
#' @param object object
#' @docType methods
#'
setMethod(f = "show", signature = "spWAVE", definition = function(object) {
  dplyr::glimpse(object)
  invisible(x = NULL)
})


#TODO 可能要改多个check_row那个函数，使得能使用多个对象
#TODO 不确定要不要做检测有无行名的前提
# #' check validity for spWAVE
# #'
# #' @param spWAVE object
# #' @param object object
# #' @docType methods
# #'
# setValidity("spWAVE", function(object) {
#   if (length(object@name) != length(object@age)) {
#     "@name and @age must be same length"
#   } else {
#     TRUE
#   }
# })


#' @inheritParams filter_LR_expr
#' @inheritParams generate_complex_data
#' @inheritParams cluster_info_identifier
create_spWAVE_object <- function(
    seurat_obj,database,
    assay="SCT",
    cluster=NULL,
    min_expr=0.1,
    min_n_cell=NULL,
    min_pct_cell=0.01,
    complex_min_cell=10
    ){
  new_object <- new("spWAVE")
  new_object@expr_raw <- 
      filter_LR_expr(db=database,seurat_obj,assay = assay,
          min_expr=min_expr,  min_n_cell=min_n_cell,   min_pct_cell=min_pct_cell
      )
  
  complex_data <- 
    generate_complex_data(database,new_object@expr_raw,
      complex_min_cell=complex_min_cell)
  
  new_object@expr_complex <- complex_data$expr_LR_df
  new_object@kept_db <- complex_data$kept_db
  new_object@coord <- get_coordinates(seurat_obj)

  new_object@cluster_info <- cluster_info_identifier(seurat_obj,cluster)
  return(new_object)
}

#*********************************
#* adapt methods for spWAVE obj **
#*********************************
#' show method for spWAVE
#'
#' @param spWAVE object
#' @param generate_kmeans_coord 
#' @param object object
#' @docType methods
#'
setMethod(f = "generate_kmeans_coord", signature = "spWAVE", definition = function(object) {
  object@meta_coord <- 
      generate_kmeans_coord(object@coord,iter.max = 15,nstart = 3)
  return(object)
})

#' show method for spWAVE
#'
#' @param spWAVE object
#' @param generate_meta_expr 
#' @param object object
#' @docType methods
#'
setMethod(f = "generate_meta_expr", signature = "spWAVE", definition = function(object) {
  object@meta_expr <- 
    generate_meta_expr(object@expr_complex,object@meta_coord$km_cluster)
  return(object)
})

#*********
#** 先对visium的方法做一个封装，并考虑一下slot的组成


setMethod(f = "perform_LR_field_calc", signature = "spWAVE", definition = function(object) {
  object@meta_expr <- 
    perform_LR_field_calc(object@kept_db,object@expr_complex,object@coord)
  return(object)
})