
#' Calc single molecular vector field
#'
#' Calc single molecular vector field from spatial coordinate and expression
#'
#' @param spatial_expr data.frame, contain spatial coordinate and expression.
#' @param gene_name selected gene name.
#' @param K_constant constant in vector field. Default is 1/(4*pi)
#'
#' @return data.frame contain, coordinates,expression, vector field
single_field_vector <- function(
    spatial_expr,
    gene_name=colnames(spatial_expr)[4],
    K_constant = 1/(4*pi)
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
