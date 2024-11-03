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
#' @export
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
#' @export
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
#' @export
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
#' @import Seurat
#' @import dplyr 
#' @return return data.frame contain sum-up vector field coordinates, expression, vector field
#' @export
perform_single_LR_spWAVE <- function(
  seurat_obj,assay="SCT",
  L_genes,R_genes
){
  #extract coordinates and expression
  calc_df <- 
    get_spatial_expr(seurat_obj,unlist(c(L_genes,R_genes)),assay = assay) %>%
    dplyr::select(c("x","y","barcode"),unlist(c(L_genes,R_genes)))
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

#*******************************
#** Interaction Score  Module **
#*******************************

#*** Spot2Spot Interaction Score

#' prepare spot2spot distance matrix
#'
#' Calc distance between every spot pairs with subtracting coordinates
#'
#' @param kept_db data.frame, LR database, must be same with the database used in db_field_result
#' @param db_field_result list, result of perform_LR_field_calc
#'
#' @return list,contains distance matrix and coordinates difference matrix of x and y
#' @export
prep_S2S_dist_mat <- function(kept_db,db_field_result){
  field_df <- concentrate_LR_field_info(
    db_field_result,
    kept_db$Ligand[1],kept_db$Receptor[1]
    )
  source_mat <- field_df[,c("x","y")] %>% as.matrix()
  pairwise_dist_mat <- 
    spa_vectorized_pdist(source_mat,source_mat)
  pairwise_sub_list <- 
    pairwised_mat_subtract(source_mat,source_mat)

  p_sub_x_mat <- do.call(cbind, lapply(pairwise_sub_list, function(mat) mat[,"x"]))
  p_sub_y_mat <- do.call(cbind, lapply(pairwise_sub_list, function(mat) mat[,"y"]))
  colnames(p_sub_x_mat) <- rownames(p_sub_x_mat)
  colnames(p_sub_y_mat) <- rownames(p_sub_y_mat)
  
  dist_mat_list <- list(p_dist=pairwise_dist_mat,p_sub_x=p_sub_x_mat,p_sub_y=p_sub_y_mat)
  return(dist_mat_list)
}


#' prepare spot2spot interaction matrix
#'
#' Prepare spot2spot interaction info including LR pair expression
#'
#' @param field_df data.frame, contain vector field info.
#' @param ligand Ligand gene.
#' @param receptor Receptor gene.
#'
#' @return return data.frame contain spot2spot interaction matrix
#' @export
prep_S2S_LR_mat <- function(
  field_df,
  ligand,receptor
){
  #** col as receiver, row as sender
  q2 <- field_df[,receptor] %>% as.matrix()
  q1 <- field_df[,ligand] %>% as.matrix()
  q_net <- q1-q2
  #q_mat <- q1 %*% t(q2)
  #q_mat <- q_net %*% t(q2)
  q_mat <- q_net %*% t(q_net)
  q_mat %<>% as("dgCMatrix")
  rownames(q_mat) <- colnames(q_mat) <- rownames(field_df)

  #** single directional checking whether the LR interaction exist
  #** 0: not exist, 1: exist
  temp <- sign(q_net) %*% t(rep(1,length(q_net)))
  LR_kept <- temp - t(temp)
  LR_kept[LR_kept<2] <- 0
  LR_kept[LR_kept==2] <- 1

  q_mat <- q_mat*LR_kept
  return(q_mat)
}


#' Prepare spot2spot interaction and distance matrix
#'
#' Prepare spot2spot interaction and distance matrix
#' A wrapper of prep_S2S_dist_mat and prep_S2S_LR_mat
#'
#' @param kept_db data.frame, LR database, must be same with the database used in db_field_result
#' @param db_field_result list, result of perform_LR_field_calc
#'
#' @return return list, the result of prep_S2S_dist_mat and prep_S2S_LR_mat
#' @export
prep_database_S2S_list <- function(kept_db,db_field_result){
  #single_mol_field_list <- db_field_result$single_mol_field_list
  #LR_pair_field_list <- db_field_result$LR_pair_field_list

  #** 在循环外进行计算距离矩阵和相互xy，后面共用

  dist_list <- prep_S2S_dist_mat(kept_db,db_field_result)

  #** for循环内每次需要提取一次LR_info,提取后再套取S2S_score进行计算
  #** 其中LR_info需要的L与R基因名也可以给到S2S_score
  LR_S2S_mat_list <- list()
  field_df_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  iter=seq_len(nrow(kept_db))
  for(i in iter){
    
    field_df_list[[i]] <- concentrate_LR_field_info(
        db_field_result,
        kept_db$Ligand[i],kept_db$Receptor[i]
        )
    
    
    LR_S2S_mat_list[[i]] <- 
      prep_S2S_LR_mat(
        field_df_list[[i]],
        kept_db$Ligand[i],kept_db$Receptor[i]
        )
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  names(LR_S2S_mat_list) <- names(field_df_list) <- kept_db$id
  prep_list <- list(dist=dist_list,q_list=LR_S2S_mat_list,field_list=field_df_list)

  return(prep_list)
}



#' calculation spot2spot interaction force
#'
#' calculation spot2spot interaction force of selected LR pair 
#' based on the distance and LR interaction matrices.
#'
#' @param dist_list distance matrix list, the result of prep_S2S_dist_mat
#' @param LR_mat_list spot2spot interaction matrix list, the result of prep_S2S_LR_mat
#' @param LR_pair selected LR pair to be calculated.
#'
#' @return list of spot2spot interaction force, including force components of x,y and norm
#' @export
calc_S2S_force_mat <- function(
  dist_list,
  LR_mat_list,
  LR_pair
){
  #** col as receiver, row as sender
  dist <- dist_list$p_dist
  sub_x <- dist_list$p_sub_x
  sub_y <- dist_list$p_sub_y
  LR_mat <- LR_mat_list[[LR_pair]]

  force_x_mat <- -LR_mat*sub_x/(dist^3)
  force_y_mat <- -LR_mat*sub_y/(dist^3)
  force_x_mat[is.na(force_x_mat)] <- 0
  force_y_mat[is.na(force_y_mat)] <- 0
  force_norm_mat <- sqrt(force_x_mat^2+force_y_mat^2)

  #force_x_mat %<>% as("dgCMatrix")
  #force_y_mat %<>% as("dgCMatrix")
  #force_norm_mat %<>% as("dgCMatrix")
  force_list <- 
    list(force_x=force_x_mat,force_y=force_y_mat,force_norm=force_norm_mat)

  return(force_list)
}


#' calculation spot2spot interaction force of database LR pairs
#'
#' calculation spot2spot interaction force of database LR pairs
#'
#' @param kept_db data.frame, LR database, must be same with the database used in prep_list
#' @param prep_list list, the result of prep_database_S2S_list
#'
#' @return list of each LR pair of spot2spot interaction force.
#' @export
calc_database_S2S_force <- function(kept_db,prep_list){
  dist_list <- prep_list$dist
  LR_mat_list <- prep_list$q_list

  db_S2S_force_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  iter=seq_len(nrow(kept_db))
  for(i in iter){
    
    db_S2S_force_list[[i]] <- 
      calc_S2S_force_mat(
        dist_list,
        LR_mat_list,
        kept_db$id[i]
        )
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  names(db_S2S_force_list) <- kept_db$id

  return(db_S2S_force_list)
}


#' calculation field-spot interaction force
#'
#' calculation field-spot interaction force which 
#' considered as how the whole field affect the spot
#'
#' @param field_df data.frame, contain vector field info.
#' @param ligand Ligand gene.
#' @param receptor Receptor gene.
#'
#' @return return list of matrices including force component of x,y and force norm.
#' @import Matrix
#' @export
calc_field_force_mat <- function(
  field_df,
  ligand,receptor
){
  #** col as receiver, row as sender
  ori_field_force <- 
  field_df %>%
    #rowwise() %>%
    mutate(
      q_net=.data[[ligand]]-.data[[receptor]],
      Fx=q_net*Ex,
      Fy=q_net*Ey
      ) %>% as.data.frame
    # mutate(
    #   Fx=.data[[receptor]]*Ex,
    #   Fy=.data[[receptor]]*Ey
    #   ) %>% as.data.frame

  #** 这里修改了一下，使用qnet进行计算而不是receptor来计算

  n_row <- rep(1,nrow(field_df)) %>% as.matrix() 
  
  field_Fx <- n_row %*% t(ori_field_force$Fx) %>% as("dgCMatrix")
  field_Fy <- n_row %*% t(ori_field_force$Fy) %>% as("dgCMatrix")
  field_Fnorm <- sqrt(field_Fx^2+field_Fy^2) %>% as("dgCMatrix")

  field_force <- 
    list(force_x=field_Fx,force_y=field_Fy,force_norm=field_Fnorm)
  field_force <- 
    lapply(field_force,function(mat){
      rownames(mat) <- colnames(mat) <- rownames(field_df)
      return(mat)
    })
  return(field_force)
}


#' calculation of spot2spot interaction score
#'
#' calculation of spot2spot interaction score of database LR pairs
#'
#' @param field_force_list field force list, the result of calc_field_force_mat.
#' @param S2S_force_list S2S force list, the result of calc_S2S_force_mat.
#' @details the interaction score is calculated by 
#' the product of components of S2S force on the field force direction and norm of field force.
#' 算了先用中文写，这一步是计算S2S force 在field force方向上的分量 乘以S2S 的模长得到的，
#' 在具体代码是线上，是直接使用了向量点乘再除以field force的模长，但两者是等价的。
#' @return return list of each LR pair of spot2spot interaction score.
#' @import Matrix
#' @export
calc_S2S_score_mat <- function(field_force_list,S2S_force_list){

  field_Fx <- field_force_list$force_x
  field_Fy <- field_force_list$force_y
  field_Fnorm <- field_force_list$force_norm
  
  S2S_Fx <- S2S_force_list$force_x
  S2S_Fy <- S2S_force_list$force_y
  S2S_Fnorm <- S2S_force_list$force_norm

  Ta <- S2S_Fx*field_Fx
  Tb <- S2S_Fy*field_Fy
  temp <- (Ta+Tb)

  S2S_cos_norm <- temp/field_Fnorm
  #S2S_cos <- temp/(field_Fnorm*S2S_Fnorm)
  #df$cos_norm <- (df$Fx*vec$Fx+df$Fy*vec$Fy)/sqrt(vec$Fx^2+vec$Fy^2)

  #** devided by 0 introducing NaN, and change the class into dgeMatrix
  S2S_cos_norm[is.na(S2S_cos_norm)] <- 0
  S2S_cos_norm %<>% as("CsparseMatrix")

  #S2S_cos[is.na(S2S_cos)] <- 0
  #S2S_cos %<>% as("CsparseMatrix")

  return(S2S_cos_norm)
}

#' calculation of spot2spot interaction score of each LR pair in database
#'
#' calculation of spot2spot interaction score of each LR pair in database,
#' a wrapper of calc_S2S_score_mat
#'
#' @param kept_db data.frame, LR database, must be same with the database used in prep_list
#' @param prep_list list, the result of prep_database_S2S_list
#' @param db_S2S_force_list list, the result of calc_database_S2S_force
#'
#' @return list of spot2spot interaction score of each LR pair in database
#' @export
calc_database_S2S_score <- function(kept_db,prep_list,db_S2S_force_list){
  dist_list <- prep_list$dist_list
  LR_mat_list <- prep_list$q_list
  field_df_list <- prep_list$field_list
  
  db_S2S_score_list <- list()
  pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  iter=seq_len(nrow(kept_db))
  for(i in iter){
    db_field_force_df <- 
      calc_field_force_mat(
        field_df_list[[i]],
        kept_db$Ligand[i],
        kept_db$Receptor[i]
        )
    
    db_S2S_score_list[[i]] <- 
      calc_S2S_score_mat(
        db_field_force_df,
        db_S2S_force_list[[i]]
        )
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  names(db_S2S_score_list) <- kept_db$id
  
  return(db_S2S_score_list)
}

#*** Cluster2Cluster Interaction Score

#' calculation of cluster2cluster interaction score matrix
#'
#' calculation of cluster2cluster interaction score matrix from sum-up S2S score by clusters
#'
#' @param S2S_score_mat S2S score matrix
#' @param clu_info_list data.frame, cluster info, including cluster name and barcode
#'
#' @return cluster2cluster interaction score matrix
#' @export
calc_C2C_mat <- function(S2S_score_mat,clu_info_list){
  n <- length(clu_info_list)
  C2C_score_mat <- matrix(nrow=n, ncol=n)
  #idx_list <- lapply(clu_info_list, function(clu) which(rownames(S2S_score_mat) %in% clu))

  for (i in 1:n) {
    tmp_idx_i <- clu_info_list[[i]]
    for (j in 1:n) {
      tmp_idx_j <- clu_info_list[[j]]
      C2C_score_mat[i,j] <- mean(S2S_score_mat[tmp_idx_i, tmp_idx_j])
    }
  }
  return(C2C_score_mat)
}

#' calculation of cluster2cluster interaction score and p-value
#'
#' calculation of cluster2cluster interaction score and p-value 
#' from sum-up S2S score by clusters and shuffle test for p-value
#'
#' @param S2S_score_mat S2S score matrix
#' @param clu_info_list data.frame, cluster info, including cluster name and barcode
#' @param clu_shuf_list data.frame, shuffled cluster info, the same format as clu_info_list.
#' Usually generated by generate_shuffle_cluster_info.
#' @return list of unfiltered cluster2cluster interaction score summary, score and p-value matrices
#' @export
calc_C2C_score_loop <- function(S2S_score_mat,clu_info_list,clu_shuf_list){
  #clac the sum force from S2S to C2C
  C2C_score <- calc_C2C_mat(S2S_score_mat,clu_info_list)

  n <- length(clu_info_list)
  C2C_p_value <- matrix(0, nrow=n, ncol=n)
  shuffle_iter <- length(clu_shuf_list)
  
  for(i in 1:shuffle_iter){
    C2C_shuf <- calc_C2C_mat(S2S_score_mat,clu_shuf_list[[i]])
    #C2C_p_value[C2C_shuf > C2C_score] <- C2C_p_value[C2C_shuf > C2C_score] + 1
    C2C_p_value <- ifelse(C2C_shuf > C2C_score, C2C_p_value + 1, C2C_p_value)
  }
  C2C_p_value <- C2C_p_value/shuffle_iter  

  rownames(C2C_score) <- names(clu_info_list)
  colnames(C2C_score) <- names(clu_info_list)
  rownames(C2C_p_value) <- names(clu_info_list)
  colnames(C2C_p_value) <- names(clu_info_list)

  score_long <- 
  C2C_score %>% 
    as.data.frame() %>%
    rownames_to_column(var="Source") %>%
    pivot_longer(!Source,names_to = "Target",values_to = "raw_score") %>%
    mutate(scale_score=scale(.data[["raw_score"]])[,1])

  p_value_long <- 
  C2C_p_value %>% 
    as.data.frame() %>%
    rownames_to_column(var="Source") %>%
    pivot_longer(!Source,names_to = "Target",values_to = "p_value")

  long_df <- full_join(score_long,p_value_long,by=c("Source","Target"))
  return(list(summary=long_df,score_mat=C2C_score,p_val_mat=C2C_p_value))
}



#' calculation of cluster2cluster interaction score of LR pairs in database
#'
#' calculation of cluster2cluster interaction score of LR pairs in database 
#' and shuffle test for p-value
#'
#' @param kept_db data.frame, LR database, must be same with the database used in db_S2S_score_list.
#' @param db_S2S_score_list database spot2spot interaction score list, the result of calc_database_S2S_score.
#' @param seurat_obj Seurat Object, extracting cluster information, optional.
#' @param cluster cluster info, see \code{\link{cluster_info_identifier}}
#' @param shuffle_iter number of shuffle iterations,default 500
#' @param random_seed random seed used in shuffle. Default is 42, 
#' the answer to the ultimate question of life, the universe, and everything.
#'
#' @return list of each LR pair C2C score summary, score and p value.
#' @export
calc_database_C2C_score <- function(
  kept_db,db_S2S_score_list,
  seurat_obj,cluster=NULL,shuffle_iter=500,random_seed=42
){
  cluster_info <- cluster_info_identifier(seurat_obj,cluster)
  #clu_info_list <- split(cluster_info$barcode,cluster_info$cluster)
  clu_info_list <- split(1:nrow(cluster_info),cluster_info$cluster)

  set.seed(random_seed)
  clu_shuf_list <- generate_shuffle_list(cluster_info,shuffle_iter)

  pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  C2C_score_list <- list()
  iter=seq_len(nrow(kept_db))
  for(i in iter){
    C2C_score_list[[i]] <- 
      calc_C2C_score_loop(
        db_S2S_score_list[[i]],
        clu_info_list,
        clu_shuf_list
      )
    setTxtProgressBar(pb, i)
  }
  close(pb)
  names(C2C_score_list) <- kept_db$id
  set.seed(NULL)

  return(C2C_score_list)
}


#' calculation of cluster2cluster interaction score of LR pairs in database
#'
#' The parallel version of calc_database_C2C_score. The parallel is based on 'future' package.
#' calculation of cluster2cluster interaction score of LR pairs in database 
#' and shuffle test for p-value
#'
#' @param kept_db data.frame, LR database, must be same with the database used in db_S2S_score_list.
#' @param db_S2S_score_list database spot2spot interaction score list, the result of calc_database_S2S_score.
#' @param seurat_obj Seurat Object, extracting cluster information, optional.
#' @param cluster cluster info, see \code{\link{cluster_info_identifier}}
#' @param shuffle_iter number of shuffle iterations,default 500
#' @param random_seed random seed used in shuffle. Default is 42, 
#' the answer to the ultimate question of life, the universe, and everything.
#'
#' @return list of each LR pair C2C score summary, score and p value.
#' @import future
#' @import future.apply
#' @import pbapply
#' @export
calc_database_C2C_score_paral <- function(
  kept_db,db_S2S_score_list,
  seurat_obj,cluster=NULL,shuffle_iter=500,random_seed=42
){
  #require(future)
  #require(future.apply)

  cluster_info <- cluster_info_identifier(seurat_obj,cluster)
  #clu_info_list <- split(cluster_info$barcode,cluster_info$cluster)
  clu_info_list <- split(1:nrow(cluster_info),cluster_info$cluster)

  set.seed(random_seed)
  clu_shuf_list <- generate_shuffle_list(cluster_info,shuffle_iter)
  #TODO 需要添加一个防呆设置检测输入s2s_force 还是S2S_score

  #handlers("progress", "beepr")
  #pb <- txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  C2C_score_list <- list()
  iter=seq_len(nrow(kept_db))
  #** future_lapply progress bar adapted from https://github.com/HenrikBengtsson/future.apply/issues/34
  #with_progress({
  #  p <- progressor(along = db_S2S_score_list)
    C2C_score_list <- 
    pblapply(db_S2S_score_list,cl="future",function(S2S_score){
        calc_C2C_score_loop(
          S2S_score,
          clu_info_list,
          clu_shuf_list
        )
      #setTxtProgressBar(pb, i)
      #p(sprintf("x=%g", S2S_score))
    })
  #})
  #close(pb)
  names(C2C_score_list) <- kept_db$id
  set.seed(NULL)

  return(C2C_score_list)
}