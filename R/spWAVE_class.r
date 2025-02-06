
#' The spWAVE Class
#'
#' The spWAVE object is created from a spatial transcriptomic Seurat V5 object.
#' When inputting an data matrix, it takes a digital data matrices as input. Genes should be in rows and cells in columns. rownames and colnames should be included.
#' The class provides functions for data preprocessing, intercellular communication network inference, communication network analysis, and visualization.
#'
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
      meta_complex = "Mat_like",
      #** field result
      single_mole_field = "list",
      LR_pair_field = "list",
      LR_family_field = "list"
      #** score result
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


#' Create spWAVE object
#' 
#' Create spWAVE object from Seurat 
#' 
#' @inheritParams filter_LR_expr
#' @inheritParams generate_complex_data
#' @inheritParams cluster_info_identifier
#' 
#' @details This function is a warper of multi helper functions, more details please check seealso link. 
#' 
#' @seealso \code{\link{cluster_info_identifier}}
#' @seealso \code{\link{generate_complex_data}}
#' @seealso \code{\link{filter_LR_expr}}
#' 
#' @return a setup spWAVE object
#' 
#' @export
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
#' Generate Meta coordinates
#'
#' Generate meta coordinates by kmeans methods
#'
#' @param coord data.frame,spatial coordinates, 
#' should have column names "x" and "y" and "barcode".
#' @param random_seed random seed of kmeans for reproducibility. Default is 42.
#' @inheritParams stats::kmeans
#' 
#' @seealso \code{\link[stats]{kmeans}}
#'
#' @return If input is data.frame, return a list, meta coordinates contain geometric center of each cluster and 
#' spot ids. And km_cluster contain the barcode and cluster id. 
#' If spWAVE, return spWAVE.
#' @export 
setGeneric("generate_kmeans_coord", function(coord,
  centers=nrow(coord)/10,
  iter.max=10,
  nstart=1,
  random_seed=42){
  standardGeneric("generate_kmeans_coord")
})


#' @rdname generate_kmeans_coord
##' @param coord object of class `spWAVE`.
#' @aliases generate_kmeans_coord,spWAVE-method
#'
#' @export
setMethod(f = "generate_kmeans_coord", signature = "spWAVE", 
  definition = function(coord,
  centers,
  iter.max=10,
  nstart=1,
  random_seed=42) {
  meta_list <- 
    generate_kmeans_coord_impl(coord=coord@coord,centers = nrow(coord@coord)/10,
      iter.max = iter.max,nstart = nstart,random_seed = random_seed)
  
  coord@meta_coord <- meta_list$meta
  coord@meta_coord_clu <- meta_list$km_cluster
  return(coord)
})

#' @rdname generate_kmeans_coord
##' @param coord object of class `data.frame`.
#' @aliases generate_kmeans_coord,data.frame-method
#' 
#' @export
setMethod("generate_kmeans_coord", "data.frame", function(coord,
  centers=nrow(coord)/10,
  iter.max=10,
  nstart=1,
  random_seed=42){
  generate_kmeans_coord_impl(
      coord,
      centers = centers,
      iter.max = iter.max,
      nstart = nstart,
      random_seed = random_seed)
})


#' Generate Meta Expression
#'
#' Generate a meta expression data.frame of clustered meta coordinates 
#' 
#' @param expr expression matrix or data.frame, commonly the result of generate_complex_data
#' @param clu_info cluster info, must be contain "barcode" and "cluster". 
#' The cluster is the result of Kmenas. 
#'
#' @details This function will sum-up the expresion of spot in a cluster. 
#' For receptor complex, we firstly calculate the complex expression of each spot by
#' generate_complex_data before this function.
#' 
#' @return data.frame, sum-uped meta expression.
#' @export
setGeneric("generate_meta_expr", function(expr,clu_info){
  standardGeneric("generate_meta_expr")
})


#' @rdname generate_meta_expr
##' @param expr object of class `spWAVE`.
#' @aliases generate_meta_expr,spWAVE-method
#'
#' @export
setMethod(f = "generate_meta_expr", signature = "spWAVE", 
  definition = function(expr,clu_info){
  expr@meta_complex <- generate_meta_expr(expr@expr_complex,expr@meta_coord_clu)
  return(expr)
})

#' @rdname generate_meta_expr
##' @param expr object of class `Mat_like`, including "matrix", "dgCMatrix","data.frame".
#' @aliases generate_meta_expr,Mat_like-method
#' 
#' @export
setMethod("generate_meta_expr", "Mat_like", function(expr,clu_info){
  generate_meta_expr_impl(expr,clu_info)
})

#' Perform calculation of LR field using hole method
#'
#' Perform calculation of LR field with hole method in 1 step.
#' 
#' @inheritParams calc_database_holed_field 
#' @inheritParams calc_database_LR_field
#' 
#' @details This function is a warper for integrating multi-steps core functions. 
#' For more details, please check the seealso.
#' 
#' @seealso \code{\link{calc_database_holed_field}}
#' @seealso \code{\link{calc_database_LR_field}}
#' 
#'
#' @export
setGeneric("perform_LR_field_hole_calc", function(
  kept_db,
  expr,
  coord,
  ROI_barcode=NULL,
  km_coord_list,
  radius = 500,
  verbose=TRUE){
  standardGeneric("perform_LR_field_hole_calc")
})

#' @rdname perform_LR_field_hole_calc
#' @aliases perform_LR_field_hole_calc,spWAVE-method
#' 
#' @export
setMethod("perform_LR_field_hole_calc", "spWAVE", function(
  kept_db,
  expr,
  coord,
  ROI_barcode=NULL,
  km_coord_list,
  radius = 500,
  verbose=TRUE){
    slot_check(kept_db,"single_mole_field")
    slot_check(kept_db,"LR_pair_field")
    slot_check(kept_db,"LR_family_field")
    meta_coord_list <- list(
      meta_coord = kept_db@meta_coord,
      km_cluster = kept_db@meta_coord_clu
    )
    print("Step1. calc single molecule or complex field")
    kept_db@single_mole_field <- calc_database_holed_field(
      kept_db=kept_db@kept_db,
      expr=kept_db@expr_complex,
      coord=kept_db@coord,
      ROI_barcode=ROI_barcode,
      km_coord_list=meta_coord_list,
      radius=radius,
      verbose=verbose
    )
    print("Step2. calc LR pair or family field")
    LR_field <- 
      calc_database_LR_field(kept_db@kept_db,kept_db@single_mole_field,verbose=verbose)
    kept_db@LR_pair_field = LR_field[[1]] 
    kept_db@LR_family_field = LR_field[[2]]
    return(kept_db)
})


#' @rdname perform_LR_field_hole_calc
#' @aliases perform_LR_field_hole_calc,Mat_like-method
#' 
#' @export
setMethod("perform_LR_field_hole_calc", "Mat_like", function(
  kept_db,
  expr,
  coord,
  ROI_barcode=NULL,
  km_coord_list,
  radius = 500,
  verbose=TRUE){
    print("Step1. calc single molecule or complex field")
    single_field <- calc_database_holed_field(
      kept_db=kept_db,
      expr=expr,
      coord=coord,
      ROI_barcode=ROI_barcode,
      km_coord_list=km_coord_list,
      radius=radius,
      verbose=verbose
    )
    print("Step2. calc LR pair or family field")
    LR_field <- calc_database_LR_field(kept_db,single_field,verbose=verbose)
    whole_db_field <- 
      list(single_mol_field_list = single_field,
          "LR_pair_field_list" = LR_field[[1]], 
          "LR_family_field_list" = LR_field[[2]])
    return(whole_db_field)
})


#' perform calculation of LR pair or family in database
#'
#' perform calculation of LR pair or family in database in 1 step.
#' 
#' @inheritParams calc_database_single_field
#' 
#' @details This function is a warper for integrating multi-steps core functions.
#' An integrated function of calc_database_single_field and calc_database_LR_field.
#' For more details, please check the seealso.
#' 
#' @seealso \code{\link{calc_database_single_field}}
#' @seealso \code{\link{calc_database_LR_field}}
#'
#' @return return list of single molecule field, 
#' LR pair field and LR family field in 3 separated list.
#' @export
setGeneric("perform_LR_field_calc", function(
  kept_db,expr,coord,verbose=TRUE){
  standardGeneric("perform_LR_field_calc")
})

#' @rdname perform_LR_field_calc
#' @aliases perform_LR_field_calc,Mat_like-method
#' 
#' @export
setMethod("perform_LR_field_calc", "Mat_like", function(
  kept_db,expr,coord,verbose=TRUE){
  print("Step1. calc single molecule or complex field")
  single_mol_field <- calc_database_single_field(kept_db,expr,coord,verbose = verbose)
  print("Step2. calc LR pair or family field")
  LR_field_list <- calc_database_LR_field(kept_db,single_mol_field,verbose = verbose)

  res_list <- list(
      single_mol_field_list=single_mol_field,
      LR_pair_field_list=LR_field_list[[1]],
      LR_family_field_list=LR_field_list[[2]]
    )
  return(res_list)
})

#' @rdname perform_LR_field_calc
#' @aliases perform_LR_field_calc,spWAVE-method
#' 
#' @export
setMethod("perform_LR_field_calc", "spWAVE", function(
    kept_db,expr,coord,verbose=TRUE){
  slot_check(kept_db,"single_mole_field")
  slot_check(kept_db,"LR_pair_field")
  slot_check(kept_db,"LR_family_field")
  print("Step1. calc single molecule or complex field")
  kept_db@single_mole_field <- 
    calc_database_single_field(kept_db@kept_db,kept_db@expr_complex,kept_db@coord,verbose = verbose)
  print("Step2. calc LR pair or family field")
  LR_field_list <- 
    calc_database_LR_field(kept_db@kept_db,kept_db@single_mole_field,verbose = verbose)

  kept_db@LR_pair_field <- LR_field_list[[1]]
  kept_db@LR_family_field <- LR_field_list[[2]]

  return(kept_db)
})

#' slot check
#'
#' check whether a slot is empty.
#'
#' @param object object of S4 class
#' @param slot slot name
#'
#' @export
slot_check <- function(object,slot){
  if(length(slot(object,slot)) !=0){
    #the slot is not empty
    warning(paste0("The slot ",slot," is not empty!"))
  }
}
