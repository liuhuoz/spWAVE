#' get_spatial_expr
#'
#' Exract coordinates and expression values from Seurat object
#'
#' @param SrtObj Seurat object.
#' @param gene one or more gene names in Seurat object.
#' @param assay assay name in slot of Seurat object. Default="RNA"
#'
#' @return return a data.frame containing coordinates, barocdes, and expression values.
#' @examples
#' # ADD_EXAMPLES_HERE
get_spatial_expr <- function(SrtObj,gene,assay="RNA"){
  coord_info <- Seurat::GetTissueCoordinates(SrtObj)
  coord_info$imagerow <- 
    max(coord_info$imagerow) - coord_info$imagerow + min(coord_info$imagerow)
  coord_info$barcode <- rownames(coord_info)

  expr_mat <- Seurat::GetAssayData(SrtObj,assay=assay)
  select_gene <- 
    expr_mat %>% as.matrix() %>% t() %>%
    as.data.frame() %>% 
    select(all_of(gene))
  colnames(select_gene) <- gene
  select_gene$barcode <- rownames(select_gene)

  spatial_expr <- dplyr::full_join(coord_info,select_gene,by="barcode")
  colnames(spatial_expr) <- c("y","x","barcode",gene)
  return(spatial_expr)
}


#' from an excellent post: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
#' this function is called by other functions to quickly compute the distance between
#' cells to grid points, or between grid points
#'
#' @param source matrix
#' @param target matrix
#' @return returns pairwise-distances
#'
spa_vectorized_pdist <- function(source,target){
  an = apply(source, 1, function(rvec) base::crossprod(rvec,rvec))
  bn = apply(target, 1, function(rvec) base::crossprod(rvec,rvec))

  m = nrow(source)
  n = nrow(target)

  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)

  #** The abs function is used to avoid the micro-tiny minus value 
  #** which produced by double precision rounding.
  #** It only occur when calc. a point to itself and will not affect others.
  mat <- abs( tmp - 2 * base::tcrossprod(source,target))
  return(sqrt(mat))
}


#' pairwise calculate 2 matrix subtraction by row between
#' cells to grid points, or between grid points,and return a list of each target row.
#' The result would be source subtracted by target.
#'
#' @param source matrix
#' @param target matrix
#' @return returns pairwise-subtract for each row list
#' 
pairwised_mat_subtract <- function(source,target){
  target_list <- lapply(seq_len(nrow(target)), function(i) target[i,,drop=FALSE])
  #expand single vector to nrow of source and preform subtract
  target_list <- lapply(target_list, function(x) rep(1,nrow(source)) %x% x)
  vec_list <- list()
  vec_list <- lapply(target_list, function(x) source-x)
  return(vec_list)
}
