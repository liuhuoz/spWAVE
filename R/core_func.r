#***************************
#** single LR calc Module **
#***************************
#' Calc single molecular vector field
#'
#' Calc single molecular vector field from spatial coordinate and expression
#'
#' @param spatial_expr data.frame, contain spatial coordinate and expression.
#' @param gene_name selected gene name.
#'
#' @return data.frame contain, coordinates,expression, vector field
single_field_vector <- function(
    spatial_expr,
    gene_name=colnames(spatial_expr)[4]
    ){
  col_idx <- which(colnames(spatial_expr)==gene_name)
  colnames(spatial_expr)[col_idx] <- "gene"

  expr_point <- spatial_expr[which(spatial_expr$gene!=0),]

  source_mat <- spatial_expr[,c("x","y"),T] %>% as.matrix()
  target_mat <- expr_point[,c("x","y"),T]  %>% as.matrix()

  spatial_dist_mat <- 
    spa_vectorized_pdist(source_mat,target_mat)
  #spatial_dist_mat[spatial_dist_mat==0] <- NA
  #min_dist <- min(spatial_dist_mat,na.rm = TRUE)
  #spatial_dist_mat[is.na(spatial_dist_mat)] <- min_dist/2
  #avoid the 0 dist which would induce Inf. 
  #1/2 min_dist thought as the spot boundary radius.
  #not so good, it seems manually add a vector to those point which may disrupt the whole vector field
  spatial_sub_list <- 
    pairwised_mat_subtract(source_mat,target_mat)

  vec_list <- 
  lapply(seq_len(length(spatial_sub_list)),function(i){
      #spatial_sub_list[[i]][spatial_sub_list[[i]]==0] <- min_dist/2
      #Also,calc the vector field at the spot boundary
      merge_df <- 
        cbind(
          spatial_sub_list[[i]],
          spatial_dist_mat[,i,drop=FALSE],
          rep(1,nrow(spatial_dist_mat)) %x% expr_point[i,col_idx,T] #only i gene in list i,and should use expr_point not all point
          )
      colnames(merge_df)[3:4] <- c("r_length","gene")
      return(merge_df)
    }
  )
# `mask old function
#  vec_list <- list()
#  for(i in 1:nrow(expr_point)){
#      vec_list[[i]] <- spatial_expr[,c("x","y"),T]-expr_point[i,c("x","y"),T]
#      vec_list[[i]]$gene <- expr_point$gene[i]
#      names(vec_list)[i] <- paste("r_vec",sep="_",i)
#      #vec_list$Qr_vec[i] <- charge_point$Q[1]*r_vec
#      #Qr_length <- sqrt(Qr_vec$x^2+Qr_vec$y^2)
#  }
# `mask old function

#note here, the x,y are not coordinate but vectors subtracted from coordinate
#dplyr::mutate efficiency isn't good, I don't know why mutate is not vectorized.
#but here I remained this part of code for better comprehension.
  vec_list <- 
  lapply(vec_list,
    function(mat){
      df <- mat %>% as.data.frame()
      df$Ex=df$gene*df$x/(df$r_length^3)
      df$Ey=df$gene*df$y/(df$r_length^3)
      df$U=df$gene/df$r_length
      df$Ex=ifelse(is.na(df$Ex),0,df$Ex) #filter out Inf which produced by the expr point itself
      df$Ey=ifelse(is.na(df$Ey),0,df$Ey)
      df$U=ifelse(is.na(df$U)|is.infinite(df$U),0,df$U)
      #df %<>% mutate(
      #  #r_length=sqrt(x^2+y^2),
      #  Ex=gene*x/(r_length^3),
      #  Ey=gene*y/(r_length^3),
      #  U=gene/r_length,
      #  Ex=ifelse(is.na(Ex),0,Ex), #filter out Inf which produced by the expr point itself
      #  Ey=ifelse(is.na(Ey),0,Ey),
      #  U=ifelse(is.na(U)|is.infinite(U),0,U)
      #)
      return(df)
    }
  )

  sigma_E_vec <- vec_list[[1]][,c("Ex","Ey","U")]
  if(length(vec_list)>1){
    for(i in 2:length(vec_list)){
      sigma_E_vec=sigma_E_vec+vec_list[[i]][,c("Ex","Ey","U")]
    }
  }

  merge_df <- cbind(spatial_expr,sigma_E_vec)
  colnames(merge_df)[col_idx] <- gene_name
  return(merge_df)
}


#' Sum-up field vector
#'
#' Sum-up field vector of selected single molecular vector of genes.
#'
#' @param vector_df_list list of data.frame, 
#' each data.frame contains a single molecular vector calculated result from `single_molecular_vector`.
#'
#' @return return data.frame contain sum-up vector field
sum_field_vector <- function(vector_df_list){
  if(is.data.frame(vector_df_list)){
      sigma_E_vec = vector_df_list[,c("Ex","Ey","U")]
    }else{
      sigma_E_vec = vector_df_list[[1]][,c("Ex","Ey","U")]
  }
  if(length(vector_df_list)>1){
    for(i in 2:length(vector_df_list)){
      sigma_E_vec=sigma_E_vec+vector_df_list[[i]][,c("Ex","Ey","U")]
    }
  }
  return(sigma_E_vec)
}


#' calculate LR pair vector
#'
#' calculate LR pair vector form single molecular field vector list
#'
#' @param gene_vec_list list of single molecular field vector
#' @param L_genes Ligand genes.
#' @param R_genes Receptor genes.
#'
#' @return return data.frame contain sum-up vector field containing coordinates, expression, vector field
calc_LR_pair_vec <- function(
  gene_vec_list,
  L_genes,R_genes
){
  #retrive used list
  L_vec_list <- gene_vec_list[L_genes]
  R_vec_list <- gene_vec_list[R_genes]

  #分割LR计算后，再计算LR对
  L_field_vec <- sum_field_vector(L_vec_list)
  R_field_vec <- sum_field_vector(R_vec_list)
  LR_field_vec <- L_field_vec - R_field_vec

  merge_spatial_df <- 
    cbind.data.frame(
      gene_vec_list[[1]][,c("x","y","barcode")],
      LR_field_vec
    )
  
  rownames(merge_spatial_df) <- merge_spatial_df$barcode
  return(merge_spatial_df)
}


#' perform calculation of single LR pair field estimation
#'
#' perform calculation of given single LR pair field estimation directly from seurat object 
#'
#' @param seurat_obj Seurat Object.
#' @param assay assay name in Seurat object. Default is "SCT".
#' @param L_genes Ligand genes.
#' @param R_genes Receptor genes.
#'
#' @return return data.frame contain sum-up vector field coordinates, expression, vector field
#' @examples
#' # ADD_EXAMPLES_HERE
perform_single_LR_spWAVE <- function(
  seurat_obj,assay="SCT",
  L_genes,R_genes
){
  #extract coordinates and expression
  calc_df <- 
    get_spatial_expr(seurat_obj,unlist(c(L_genes,R_genes)),assay = assay) %>%
    select(c("x","y","barcode"),unlist(c(L_genes,R_genes)))
  #generate gene expression list and divide into Ligand and Receptor.
  gene_df_list <- 
  lapply(calc_df[,-c(1:3)],
    function(gene) cbind.data.frame(spatial_expr[,c("x","y","barcode")],gene)
  )
  for(i in seq_len(length(gene_df_list))){
    colnames(gene_df_list[[i]])[4] <- names(gene_df_list)[i]
  }
  #calc single molecular vector
  gene_vec_list <- 
    lapply(gene_df_list, function(x) single_field_vector(x))
  names(gene_vec_list) <- names(gene_df_list)
  #calc LR pair vector
  merge_spatial_df <- calc_LR_pair_vec(gene_vec_list,L_genes,R_genes)
  return(merge_spatial_df)
}


#' calculate field strength and other quantities
#'
#' calculate field strength and other quantities based on the vector field
#'
#' @param vector_df data.frame, 
#' vector field indicated by "Ex","Ey", and optionally "U".
#' @param K_constant numeric, the constant of field transmission, default is 1/(4*pi).
#'
#' @return extended data.frame of vector_df, added field strength. 
calc_field_strength <- function(vector_df,K_constant = 1/(4*pi)){
    vector_df %<>%
      mutate(
        KEx=K_constant*Ex,
        KEy=K_constant*Ey,
        E_strength=ifelse(is.infinite(Ex*Ey),Inf,K_constant*sqrt(Ex^2+Ey^2)),
      )
    if("U" %in% colnames(vector_df)){
      vector_df %<>%
        mutate(KU=K_constant*U)
    }
    return(vector_df)
}

#*****************************
#** database LR calc Module **
#*****************************
#** 目前database 计算模块返回的是情况仍然是各类零散的list，只是为了下一步引入S4对象做准备
#** 因此暂时不会添加perform_LR_database_spWAVE的整合函数


#' calculate single molecule field in database
#'
#' calculate single molecule fields of all genes in database
#'
#' @param kept_db  data.frame, LR database.
#' @param expr expression matrix.
#' @param coord spatial coordinates. 
#'
#' @return list of single molecule field
#' @examples
#' #prepare database and expression data
#' mouse_db <-
#'   load_database(db_source = "CellChat",db_species = "mouse",filter_type = "Secreted")
#' expr_ft_mat <-
#'   filter_LR_expr(db=mouse_db,brain,assay = "SCT")
#' complex_data <- generate_complex_data(mouse_db,expr_ft_mat)
#' #extract coordinates from seurat object
#' brain_coord <- get_coordinates(brain)
#' # perform calculation
#' brain_single_mol_field <- 
#'  calc_database_single_field(
#'    kept_db = complex_data$kept_db,
#'    expr = complex_data$expr_LR_df,
#'    coord = brain_coord)
calc_database_single_field <- function(kept_db,expr,coord){ 
  single_mol_field <- list()
  #print("Step1. calc single molecule or complex field")
  pb <- txtProgressBar(min = 0, max = ncol(expr), style = 3)
  for(i in seq_len(ncol(expr))){
    spatial_expr <- 
      concentrate_coord_expr(expr=expr,genes=colnames(expr)[i],coord_info=coord)
    single_mol_field[[i]] <- 
      single_field_vector(spatial_expr,gene_name=colnames(expr)[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(single_mol_field) <- colnames(expr)
  return(single_mol_field)
}

#' calculate LR pair field in database
#'
#' calculate LR pair fields by sum-up single molecular field in list and database
#'
#' @param kept_db data.frame, LR database, must be same with the database used in single_mol_field
#' @param single_mol_field list of single molecular field, 
#' usually the result of calc_database_single_field 
#'
#' @return return list of LR pair or family field estimation.
#' @examples
#' #upstream step please see calc_database_single_field.
#' brain_db_field <- 
#'   calc_database_LR_field(complex_data$kept_db,brain_single_mol_field)
calc_database_LR_field <- function(kept_db,single_mol_field){
  LR_pair_field_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  for(i in seq_len(nrow(kept_db))){
    LR_pair_field_list[[i]] <- 
      calc_LR_pair_vec(single_mol_field,kept_db$Ligand[i],kept_db$Receptor[i])
    LR_pair_field_list[[i]]$LR_pair <- 
      paste(kept_db$Ligand[i],kept_db$Receptor[i],sep=".")
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(LR_pair_field_list) <- kept_db$id

  family_lig_list <- split(kept_db$Ligand, kept_db$Family)
  family_rec_list <- split(kept_db$Receptor, kept_db$Family)
  names(family_lig_list) <- 
    names(family_rec_list) <- 
      kept_db$Family %>% unique %>% sort

  family_lig_list %<>% lapply(unique)
  family_rec_list %<>% lapply(unique)

  LR_family_field_list <- list()
  pb <- txtProgressBar(min = 0, max = length(family_lig_list), style = 3)
  for(i in seq_len(length(family_lig_list))){
    LR_family_field_list[[i]] <- 
      calc_LR_pair_vec(single_mol_field,family_lig_list[[i]],family_rec_list[[i]])
    LR_family_field_list[[i]]$Family <- names(family_lig_list)[i]
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(LR_family_field_list) <- names(family_lig_list)
  
  LR_field_list <- 
    list(LR_pair_field_list=LR_pair_field_list,
        LR_family_field_list=LR_family_field_list)
  return(LR_field_list)
}


#' perform calculation of LR pair or family in database
#'
#' perform calculation of LR pair or family in database in 1 step.
#' An integrated function of calc_database_single_field and calc_database_LR_field
#'
#' @param kept_db  data.frame, LR database.
#' @param expr expression matrix.
#' @param coord spatial coordinates. 
#'
#' @return return list of single molecule field, 
#' LR pair field and LR family field in 3 separated list.
#' @examples
#' #prepare database and expression data
#' mouse_db <-
#'   load_database(db_source = "CellChat",db_species = "mouse",filter_type = "Secreted")
#' expr_ft_mat <-
#'   filter_LR_expr(db=mouse_db,brain,assay = "SCT")
#' complex_data <- generate_complex_data(mouse_db,expr_ft_mat)
#' #extract coordinates from seurat object
#' brain_coord <- get_coordinates(brain)
#' # perform calculation
#' brain_db_res <- 
#'  perform_LR_field_calc(
#'    kept_db = complex_data$kept_db,
#'    expr = complex_data$expr_LR_df,
#'    coord = brain_coord)
perform_LR_field_calc <- function(kept_db,expr,coord){
  print("Step1. calc single molecule or complex field")
  single_mol_field <- calc_database_single_field(kept_db,expr,coord)
  print("Step2. calc LR pair or family field")
  LR_field_list <- calc_database_LR_field(kept_db,single_mol_field)

  res_list <- list(
      single_mol_field_list=single_mol_field,
      LR_pair_field_list=LR_field_list[[1]],
      LR_family_field_list=LR_field_list[[2]]
    )
  return(res_list)
}