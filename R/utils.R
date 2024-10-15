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

#' spa_vectorized_pdist
#' 
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


#' pairwised_mat_subtract
#' 
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



#' Load Ligand and Receptor Database
#'
#' load database from built-in data, and filter specific type LR pairs.
#'
#' @param db_source database source, including "CellChat" and "CellPhoneDB".
#' @param db_species species,including "human", "mouse", and "zebrafish".
#' @param filter_type kept LR type in result, including "Contact","ECM","Secreted". Default is "None", and all would be kept.
#'
#' @return return a data.frame after filtered
#' @examples
#' # cpdb_human <- load_database(db_source="CellPhoneDB",db_species="human")
load_database <- function(
  db_source=c("CellChat","CellPhoneDB"),
  db_species=c("human","mouse","zebrafish"),
  filter_type="None"
){
  if(db_source=="CellPhoneDB" & db_species=="zebrafish"){
    warning(
      "CellPhoneDB does not contain database of zebrafish, using CellChat database instead."
      )
    db_source <- "CellChat"
  }
  db_source <- match.arg(db_source)
  db_species <- match.arg(db_species)
  db_source <- ifelse(db_source=="CellPhoneDB","cpdb","cellchat")

  db_attr <- paste(db_source,db_species,sep="_")
  db <- get(db_attr)

  if(!(filter_type %in% c("Contact","ECM","Secreted","None"))){
    warning("filter_type not found in crurated type, using default 'None' instead.")
    filter_type <- "None"
  }

  if(filter_type=="None"){
      filter_type <- c("Contact","ECM","Secreted")
    }
  #** filter out same ligand and receptor which is not considered
  temp <- which(db$Ligand==db$Receptor)
  if(length(temp)>0){db <- db[-temp,]}
  db_ft <- db %>%
    filter(Type2 %in% c(filter_type))
  return(db_ft)
}


#' Filter Ligand and Receptor expression 
#' 
#' Extract expression matrix from Seurat object, and filtered
#' by ligand and receptor genes in database with expression level.
#'
#' @param db database data.frame.
#' @param SrtObj Seurat object.
#' @param assay assay slot in Seurat object. default "RNA".
#' @param min_expr minimum expression. default 0.1
#' @param min_n_cell minimum number of cell expressing a gene. default NULL. This arg will mask 'min_pct_cell'.
#' @param min_pct_cell minimum percentage of cell expressing a gene. 
#'
#' @return filtered expression matrix

filter_LR_expr <- function(
  db,
  SrtObj,
  assay="RNA",
  min_expr=0.1,
  min_n_cell=NULL,
  min_pct_cell=0.01
){
  expr_mat <- GetAssayData(SrtObj,assay=assay)
  db_all_genes <- 
    c(db$Ligand,str_split(db$Receptor,pattern = "_")) %>% 
    unlist() %>% unique()

  #filterd by database
  intersect_genes <- intersect(db_all_genes,rownames(expr_mat))
  expr_mat <- expr_mat[intersect_genes,]
  
  #fitlerd by expr
  expr_mat_ft <- filter_expr_by_cutoff(expr_mat,
    min_expr=min_expr,
    min_n_cell=min_n_cell,
    min_pct_cell=min_pct_cell
  )
  return(expr_mat_ft)
}

#' filter_expr_by_cutoff
#' filter expression matrix by expression level
#'
#' @param expr_mat expression matrix.
#' @param min_expr minimum expression.
#' @param min_n_cell minimum number of cell expressing a gene. This arg will mask 'min_pct_cell'.
#' @param min_pct_cell minimum percentage of cell expressing a gene. 
#' @return filtered expression matrix.

filter_expr_by_cutoff <- function(
  expr_mat,min_expr=0.1,
  min_n_cell=NULL,
  min_pct_cell=0.01
){
  min_n_cell <- ifelse(is.null(min_n_cell),min_pct_cell*ncol(expr_mat),min_n_cell)
  expr_mat_binary <- expr_mat %>% as.matrix()
  expr_mat_binary[expr_mat_binary < min_expr] <- 0
  expr_mat_binary[expr_mat_binary >= min_expr] <- 1
  kept_genes <- rowSums(expr_mat_binary) >= min_n_cell
  expr_mat_ft <- expr_mat[kept_genes,]
  return(expr_mat_ft)
}


#' Generate Receptor Complex Expression 
#'
#' Generate receptor complex expression by expression of each sub-unit of the complex. 
#' And aggregate the expression and database result for next step.
#'
#' @param db LR database dataframe
#' @param expr_mat gene expression matrix
#' @param complex_min_cell minimum number of cell expressing complex
#' 
#' @details The receptor complex expression are calculated using 
#' \deqn{R_{complex} = \prod_{i \in R} R_i ^{1/n_{R}}} where 
#' R is the expression of each sub-unit and n_{R} is the number of sub-units in
#' the complex.
#' 
#' @return Ruturn a list including filtered database, LR expression list, and merged LR expression dataframe.

generate_complex_data <- function(db,expr_mat,complex_min_cell=10){
  expr_df <- expr_mat %>% as.matrix() %>% t() %>% as.data.frame()
  expr_LR_list <- list()
  kept_db_index <- c()
  for(i in seq_len(nrow(db))){
    lig <- db$Ligand[i]
    rec_complex <- db$Receptor[i]
    rec <- str_split_1(rec_complex,pattern = "_")
    if(all(c(lig,rec) %in% colnames(expr_df))){ #filter out not-exist in expr LR_pair
      LR_complex_df <- expr_df[,c(lig,rec)]
      temp <- rowProds(LR_complex_df[,c(rec)] %>% as.matrix())^(1/length(rec))
      if(length(which(temp!=0))>=complex_min_cell){ #filter out low expr LR_pair
        LR_complex_df$complex <- temp
        colnames(LR_complex_df)[length(c(lig,rec,1))] <- rec_complex
        expr_LR_list[[i]] <- LR_complex_df
        kept_db_index[i] <- i
      }else{
        next
      }
    }else{
      next
    }
  }
  expr_LR_list <- expr_LR_list[!sapply(expr_LR_list,is.null)]

  temp_df <- do.call(cbind, expr_LR_list)
  dup_cols <- duplicated(colnames(temp_df))
  expr_LR_df <- temp_df[, !dup_cols]

  kept_db_index <- kept_db_index[!is.na(kept_db_index)]
  kept_db <- db[kept_db_index,]
  kept_db$id <- paste(kept_db$Ligand,kept_db$Receptor,sep=".")
  return(list(kept_db=kept_db,expr_LR_list=expr_LR_list,expr_LR_df=expr_LR_df))
}