#***************************
#** single LR calc Module **
#***************************
#' Calc single molecular vector field
#'
#' Calc single molecular vector field from spatial coordinate and expression
#'
#' @param spatial_expr data.frame, contain spatial coordinate and expression.
#' @param gene_name selected gene name.
#' @param dist_mat pre-calced pairwise distance matrix. Default is NULL.
#' If NULL, it will be calculated automatically.
#' we recommend to pre-calculate distance matrix for speed up.
#' @param coord_sub_list pre-calced pairwise subtraction matrix. Default is NULL.
#' Treatment and recommendation likes parameter \code{dist_mat}.
#'
#' @return data.frame contain, coordinates,expression, vector field
#' @export
single_field_vector <- function(
    spatial_expr,
    gene_name=colnames(spatial_expr)[4],
    dist_mat=NULL,
    coord_sub_list=NULL
    ){
  col_idx <- which(colnames(spatial_expr)==gene_name)
  colnames(spatial_expr)[col_idx] <- "gene"
  idx_temp <- which(spatial_expr$gene!=0)
  expr_point <- spatial_expr[idx_temp,]

  source_mat <- spatial_expr[,c("x","y"),T] %>% as.matrix()
  target_mat <- expr_point[,c("x","y"),T]  %>% as.matrix()
  
  if(is.null(dist_mat)){
    spatial_dist_mat <- 
      spa_vectorized_pdist(source_mat,target_mat)
  }else{
    spatial_dist_mat <- dist_mat[,idx_temp,drop=FALSE]
  }
  #spatial_dist_mat[spatial_dist_mat==0] <- NA
  #min_dist <- min(spatial_dist_mat,na.rm = TRUE)
  #spatial_dist_mat[is.na(spatial_dist_mat)] <- min_dist/2
  #avoid the 0 dist which would induce Inf. 
  #1/2 min_dist thought as the spot boundary radius.
  #not so good, it seems manually add a vector to those point which may disrupt the whole vector field
  if(is.null(coord_sub_list)){
    spatial_sub_list <- 
      pairwised_mat_subtract(source_mat,target_mat)
  }else{
    spatial_sub_list <- coord_sub_list[idx_temp]
  }
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

  merge_spatial_df <- LR_field_vec
    # cbind.data.frame(
    #   gene_vec_list[[1]][,c("x","y","barcode")],
    #   LR_field_vec
    # )
  
  rownames(merge_spatial_df) <- rownames(merge_spatial_df)
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
    function(gene) cbind.data.frame(calc_df[,c("x","y","barcode")],gene)
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
#' @param skip_subunit logical, default is TURE and will not calc fields of receptor subunit.
#' @inheritParams auto_select_lapply
#'
#' @return list of single molecule field
#' @export
calc_database_single_field <- function(kept_db,expr,coord,skip_subunit=TRUE,verbose=TRUE){ 
  single_mol_field <- list()
  if(skip_subunit){
    expr <- filter_subunit(kept_db=kept_db,expr=expr)
  }
  spatial_expr_list <- 
    lapply(colnames(expr),function(gene_name){
      concatenate_coord_expr(expr=expr,genes=gene_name,coord_info=coord)
      }
    )
  coord_mat <- coord[,c("x","y")] %>% as.matrix()
  dist_mat <- spa_vectorized_pdist(coord_mat,coord_mat)
  coord_sub_list <- pairwised_mat_subtract(coord_mat,coord_mat)
  single_mol_field <- 
    auto_select_lapply(spatial_expr_list,function(spatial_expr){
      single_field_vector(spatial_expr = spatial_expr,
      dist_mat = dist_mat,
      coord_sub_list = coord_sub_list
      )
    },verbose=verbose
    )
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
#' @inheritParams auto_select_lapply
#'
#' @return return list of LR pair or family field estimation.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
calc_database_LR_field <- function(kept_db,single_mol_field,verbose=TRUE){
  LR_pair_field_list <- list()
  if(verbose){
    pb <- utils::txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  }
  for(i in seq_len(nrow(kept_db))){
    LR_pair_field_list[[i]] <- 
      calc_LR_pair_vec(single_mol_field,kept_db$Ligand[i],kept_db$Receptor[i])
    # LR_pair_field_list[[i]]$LR_pair <- 
    #   paste(kept_db$Ligand[i],kept_db$Receptor[i],sep=".")
    if(verbose){utils::setTxtProgressBar(pb, i)}
  }
  if(verbose){close(pb)}
  names(LR_pair_field_list) <- kept_db$id

  family_lig_list <- split(kept_db$Ligand, kept_db$Family)
  family_rec_list <- split(kept_db$Receptor, kept_db$Family)
  names(family_lig_list) <- 
    names(family_rec_list) <- 
      kept_db$Family %>% unique %>% sort

  family_lig_list %<>% lapply(unique)
  family_rec_list %<>% lapply(unique)

  LR_family_field_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = length(family_lig_list), style = 3)
  for(i in seq_len(length(family_lig_list))){
    LR_family_field_list[[i]] <- 
      calc_LR_pair_vec(single_mol_field,family_lig_list[[i]],family_rec_list[[i]])
    #LR_family_field_list[[i]]$Family <- names(family_lig_list)[i]
    utils::setTxtProgressBar(pb, i)
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
#' @inheritParams calc_database_single_field 
#'
#' @return return list of single molecule field, 
#' LR pair field and LR family field in 3 separated list.
#' @export
perform_LR_field_calc_impl <- function(kept_db,expr,coord,verbose=TRUE){
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
}

#**** Centroid calc module *****
#' Calculate a single point vector
#'
#' Calculate a single point vector from spatial coordinates, 
#' for the centroid method ST like 10x Xenium 
#'
#' @param X,Y The coordinates of the spot. 
#' @param spatial_expr the spatial expression 
#' containing spatial coordinates barcode, and expression. 
#' @param dist_mat the pre-calced distance matrix, default is NULL.
#' And it will be calculated automatically if NULL. 
#' For large matrix, we recommend to auto generate the distance matrix in this function.
#' Because generating one-column matrix is faster than slicing the large fullback matrix.
#' @param coord_sub_list the pre-calced pairwise subtraction matrix, default is NULL.
#' Treatment and recommendation likes parameter \code{dist_mat}.
#'
#' @details This function is blahblah
#' 
#' 
#' @return data.frame, each row is a single point field vector of a gene 
#' at given coordinates.
#' @export 
single_point_vector <- function(
    X,Y,
    spatial_expr,
    #gene_name=colnames(spatial_expr)[4],
    dist_mat=NULL,
    coord_sub_list=NULL
  ){
  #col_idx <- which(colnames(spatial_expr)==gene_name)
  #colnames(spatial_expr)[col_idx] <- "gene"

  expr_point <- spatial_expr#[which(spatial_expr$gene!=0),]

  source_mat <- matrix(data=c(X,Y),nrow=1)
  colnames(source_mat) <- c("x","y")
  target_mat <- expr_point[,c("x","y"),TRUE] %>% as.matrix()

  if(is.null(dist_mat)){
    spatial_dist_mat <- 
      spa_vectorized_pdist(target_mat,source_mat) 
    #** here is one point to others,for returning one column,
    #** I reverse the source and target in args and return as a transposed the result.
  }else{
    spatial_dist_mat <- dist_mat
  }

  if(is.null(coord_sub_list)){
    spatial_sub_mat <- 
      rep(1,nrow(target_mat)) %x% source_mat - target_mat
  }else{
    spatial_sub_mat <- coord_sub_list
  }
  
  expr_point <- expr_point[,-c(1:3)]
  genes <- colnames(expr_point)
  vec <- 
  cbind(spatial_sub_mat,
        spatial_dist_mat,
        expr_point
    )
  colnames(vec)[3] <- c("r_length")

#note here, the x,y are not coordinate but vectors subtracted from coordinate
  vec %<>% as.data.frame() 
  vec$Ex=vec$x/(vec$r_length^3)
  vec$Ey=vec$y/(vec$r_length^3)
  vec$Ex=ifelse(is.na(vec$Ex),0,vec$Ex) #filter out Inf which produced by the expr point itself
  vec$Ey=ifelse(is.na(vec$Ey),0,vec$Ey)

  E_vec_list <- 
  lapply(genes,function(one_gene){
    E_vec <- vec
    E_vec$Ex=E_vec[[one_gene]]*E_vec$Ex
    E_vec$Ey=E_vec[[one_gene]]*E_vec$Ey
    E_vec$U=E_vec[[one_gene]]/vec$r_length
    E_vec$U=ifelse(is.na(E_vec$U)|is.infinite(E_vec$U),0,E_vec$U)
    sigma_E_vec <- 
      E_vec[,c("Ex","Ey","U")] %>% 
      colSums %>% 
      matrix(nrow=1,dimnames = list(NULL,c("Ex","Ey","U")))
      merge_df <- cbind(source_mat,sigma_E_vec) %>% as.data.frame()
  })
  gene_point_vector <- do.call(rbind,E_vec_list)
  rownames(gene_point_vector) <- genes
  #colnames(merge_df)[col_idx] <- gene_name
  return(gene_point_vector)
}

#' Calculate a single point vector
#'
#' Calculate a single point vector from spatial coordinates, 
#' for the centroid method ST like 10x Xenium 
#' 
#' @inheritParams calc_database_single_field
#' @inheritParams generate_holed_coord_expr
#' @inheritParams auto_select_lapply
#' @param ROI_barcode character vector, the barcode of spots/cells in the region of interest.
#' @param skip_subunit logical, default is TURE and will not calc fields of receptor subunit.
#'
#' @details This function is wrapper of \code{\link{single_point_vector}} 
#' and \code{\link{generate_holed_coord_expr}} to generate a database field vector of 
#' each single spot or point in spatial by using the holed methods.
#' Although the inputs is similar to \code{\link{calc_database_single_field}},
#' this function uses different method and automatically processed the data  
#' during the calculation.
#' 
#' @seealso \code{\link{generate_holed_coord}}
#' @seealso \code{\link{single_point_vector}}
#' 
#' @return list of data.frame, 
#' each data.frame is a single molecule field vector of a gene
#' @export 
calc_database_holed_field <- function(
  kept_db,
  expr,
  coord,
  ROI_barcode=NULL,
  km_coord_list,
  radius = 500,
  skip_subunit=TRUE,
  verbose=TRUE
){
  #* 这里就是直接对每个i点直接生成merge的空间表达后，直接计算，求和
  #meta_coord <- km_coord_list[[1]]
  spot_coord <- coord
  
  spot_expr <- expr
  #spot_expr$cluster <- rownames(spot_expr)
  meta_expr <- generate_meta_expr(expr,km_coord_list$km_cluster)
  all_expr <- rbind.data.frame(spot_expr,meta_expr)
  if(skip_subunit){
    all_expr <- filter_subunit(kept_db=kept_db,expr=all_expr)
  }
  if(is.null(ROI_barcode)){ROI_barcode <- rownames(spot_expr)}

  point_vec_list <- 
  auto_select_lapply(seq_len(length(ROI_barcode)),
    function(i){
    temp <- generate_holed_coord_expr(
      coord=coord,
      km_coord_list = km_coord_list,
      all_expr=all_expr,
      center=ROI_barcode[i],
      radius=radius)
    center_idx <- which(temp$barcode == ROI_barcode[i])
    temp_mat <- temp[,-3] %>% as.matrix()
    nonzero_expr <- which(rowSums(temp_mat[,3:ncol(temp_mat)])!=0)
    point_vec <- 
      single_point_vector_rcpp(temp$x[center_idx],temp$y[center_idx],temp_mat[nonzero_expr,])
    point_vec$barcode <- ROI_barcode[i]
    point_vec$x <- temp$x[center_idx]
    point_vec$y <- temp$y[center_idx]
    return(point_vec)
  },verbose=verbose
  )
  spot_expr_ROI <- spot_expr[ROI_barcode,]
  char_index <- point_vec_list[[1]]$gene
  gene_vec_array <- array(
      dim = c(length(ROI_barcode), 4, length(char_index)),
      dimnames = list(ROI_barcode, c("expression", "Ex", "Ey", "U"), char_index)
    )

  field_mat <- do.call(rbind, lapply(point_vec_list, function(df) {
    df[1, c("Ex", "Ey", "U")]  # 取第一行因为所有行坐标相同
  }))
  rownames(field_mat) <- sapply(point_vec_list, function(df) df$barcode[1])

  for(i in seq_along(char_index)) {
    gene_vec_array[, "expression", i] <- spot_expr_ROI[ROI_barcode, char_index[i]]
    gene_vec_array[,  c("Ex", "Ey", "U"), i] <- as.matrix(field_mat)
  }

  return(gene_vec_array)
}

#******************************
#** Interaction Score Module **
#******************************

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
  field_df <- concatenate_LR_field_info(
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
  
  dist_mat_list <- list(p_dist=pairwise_dist_mat,p_sub_x=p_sub_x_mat, p_sub_y=p_sub_y_mat)
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
#' @importFrom methods as
#' @importFrom MatrixExtra t_shallow
#' @importClassesFrom Matrix dgCMatrix
#' @export
prep_S2S_LR_mat <- function(
  field_df,
  ligand,receptor
){
  #** col as receiver, row as sender
  q2 <- field_df[,receptor]
  q1 <- field_df[,ligand]
  q_net <- (q1-q2) %>% as.matrix() %>% methods::as("dgCMatrix")
  #q_mat <- q1 %*% t(q2)
  #q_mat <- q_net %*% t(q2)
  #print(class(q_net))
  #return(q_net)
  q_mat <- q_net %*% MatrixExtra::t_shallow(q_net)
  #q_mat %<>% as("dgCMatrix")
  rownames(q_mat) <- colnames(q_mat) <- rownames(field_df)

  #** single directional checking whether the L to R interaction exist
  #** col to row ,col as sender, row as receiver
  #** 0: not exist, 1: exist
  #print("prep sign")
  q_net_sign <- q_net
  q_net_sign@x <- sign(q_net_sign@x)
  LR_kept <- compute_LR_kept_dgC(q_net_sign)

  q_mat <- q_mat * LR_kept #** 这里到也可以用@x来计算，不过得先对q_mat先进行处理过滤ij数值
  return(q_mat)
}

#' prepare spot2spot sparse distance matrix
#'
#' Calc distance between selected spot and return sparse matrix
#'
#' @param LR_mat LR interaction matrix, the result of prep_S2S_LR_mat
#' @param coord data.frame, coordinates of spots, 
#' the rownames should be same as rownames of LR_mat
#'
#' @return dgCMatrix, the distance between selected spot
#' @export
prep_S2S_dist_sparse <- function(LR_mat,coord){
  #align row index between LR_mat and coord
  if(all(rownames(coord) %in% rownames(LR_mat))){
    if(!all(rownames(coord) == rownames(LR_mat))){
      coord <- coord[rownames(LR_mat),]
    }
  }else(
    stop("rownames of coord should be same as rownames of LR_mat")
  )
  coord <- coord[,c("x","y")]

  LR_mat %<>% methods::as("TsparseMatrix")
  row_idx <- LR_mat@i
  col_idx <- LR_mat@j

  row_point <- coord[row_idx+1,]
  col_point <- coord[col_idx+1,]
  dist_sub <- (row_point-col_point)
  dist <- sqrt(dist_sub^2 %>% rowSums())

  dist_mat <- dist_sub_x <- dist_sub_y <- LR_mat
  dist_mat@x <- dist
  dist_sub_x@x <- dist_sub[,c("x")]
  dist_sub_y@x <- dist_sub[,c("y")]
  #dist_mat %<>% methods::as("CsparseMatrix")
  #dist_sub_x %<>% methods::as("CsparseMatrix")
  #dist_sub_y %<>% methods::as("CsparseMatrix") 
  #** 这个地方drop0了，导致长度不一致，（也可能不需要转换），想想怎么处理
  #** 还有一种就是手动算ij了，现在先不进行转换也不进行drop0，看看怎么样 

  return(
    list(p_dist=dist_mat,
    p_sub_x=dist_sub_x,
    p_sub_y=dist_sub_y)
  )
}

#' Prepare spot2spot interaction and distance matrix
#'
#' Prepare spot2spot interaction and distance matrix
#' A wrapper of prep_S2S_dist_mat and prep_S2S_LR_mat
#'
#' @param kept_db data.frame, LR database, must be same with the database used in db_field_result
#' @param db_field_result list, result of perform_LR_field_calc
#' @param coord data.frame, coordinates of spots, default is null and extract form db_field_result
#'
#' @return return list, the result of prep_S2S_dist_mat and prep_S2S_LR_mat
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
prep_database_S2S_list <- function(kept_db,db_field_result,coord=NULL){
  #single_mol_field_list <- db_field_result$single_mol_field_list
  #LR_pair_field_list <- db_field_result$LR_pair_field_list

  #** 在循环外进行计算距离矩阵和相互xy，后面共用

  #dist_list <- prep_S2S_dist_mat(kept_db,db_field_result)
  if(is.null(coord)){
    coord <- db_field_result[[1]][[1]][,c("x","y")]
  }
  #** for循环内每次需要提取一次LR_info,提取后再套取S2S_score进行计算
  #** 其中LR_info需要的L与R基因名也可以给到S2S_score
  LR_S2S_mat_list <- list()
  LR_dist_list <- list()
  field_df_list <- list()

  pb <- utils::txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  iter=seq_len(nrow(kept_db))
  # dist_flag <- TRUE
  # if(nrow(coord)<=13000){
  #   dist_list <- prep_S2S_dist_mat(kept_db,db_field_result)
  #   LR_dist_list <- rep(list(dist_list),nrow(kept_db))
  #   dist_flag <- FALSE
  # } #** 为了后面直接使用稀疏矩阵计算，这里全部使用后续的计算方式
  for(i in iter){
    
    field_df_list[[i]] <- concatenate_LR_field_info(
        db_field_result,
        kept_db$Ligand[i],kept_db$Receptor[i]
        )
    
    LR_S2S_mat_list[[i]] <- 
      prep_S2S_LR_mat(
        field_df_list[[i]],
        kept_db$Ligand[i],kept_db$Receptor[i]
        )
        
    #if(dist_flag){
      LR_dist_list[[i]] <- 
        prep_S2S_dist_sparse(
          LR_S2S_mat_list[[i]],
          coord
          )
    #}

    utils::setTxtProgressBar(pb, i)
    
  }
  close(pb)
  names(LR_S2S_mat_list) <- names(field_df_list) <- kept_db$id
  prep_list <- 
    list(dist=LR_dist_list,
      q_list=LR_S2S_mat_list,
      field_list=field_df_list)

  return(prep_list)
}



#' calculation spot2spot interaction force
#'
#' calculation spot2spot interaction force of selected LR pair 
#' based on the distance and LR interaction matrices.
#'
#' @param dist_list distance matrix list, the result of prep_S2S_dist_mat
#' @param LR_mat spot2spot interaction matrix of selected LR pair, in the result of prep_S2S_LR_mat
#' 
#' @importFrom Matrix drop0
#' 
#' @return list of spot2spot interaction force, including force components of x,y and norm
#' @export
calc_S2S_force_mat <- function(
  dist_list,
  #LR_mat_list,
  LR_mat
){
  #** col as receiver, row as sender
  dist <- dist_list$p_dist
  sub_x <- dist_list$p_sub_x
  sub_y <- dist_list$p_sub_y
  #LR_mat <- LR_mat_list[[LR_pair]]


  force_norm_mat <- force_x_mat <- force_y_mat <- LR_mat

  force_x_mat@x <- -LR_mat@x*sub_x@x/(dist@x^3)
  force_y_mat@x <- -LR_mat@x*sub_y@x/(dist@x^3)

  force_x_mat@x[is.na(force_x_mat@x)] <- 0
  force_y_mat@x[is.na(force_y_mat@x)] <- 0
  force_norm_mat@x <- sqrt(force_x_mat@x^2+force_y_mat@x^2)

  #force_x_mat %<>% Matrix::drop0()
  #force_y_mat %<>% Matrix::drop0()
  #force_norm_mat %<>% Matrix::drop0()

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
#' @inheritParams auto_select_lapply
#'
#' @return list of each LR pair of spot2spot interaction force.
#' @export
calc_database_S2S_force <- function(kept_db,prep_list,verbose=TRUE){
  dist_list <- prep_list$dist
  LR_mat_list <- prep_list$q_list

  db_S2S_force_list <- list()

    db_S2S_force_list <- 
      auto_select_lapply(seq_len(nrow(kept_db)),function(i){
        calc_S2S_force_mat(
          dist_list[[i]],
          LR_mat_list[[i]])
        },verbose=verbose)
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
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom methods as
#' @export
calc_field_force_mat <- function(
  field_df,
  ligand,receptor
){
  #** col as receiver, row as sender
  # ori_field_force <- 
  # field_df %>%
  #   #rowwise() %>%
  #   mutate(
  #     q_net=.data[[ligand]]-.data[[receptor]],
  #     Fx=q_net*Ex,
  #     Fy=q_net*Ey
  #     ) %>% as.data.frame
    field_df$q_net <- field_df[,ligand] - field_df[,receptor]
    field_df$Fx=field_df$q_net*field_df$Ex
    field_df$Fy=field_df$q_net*field_df$Ey
    ori_field_force <- field_df %>% as.data.frame()
    # mutate(
    #   Fx=.data[[receptor]]*Ex,
    #   Fy=.data[[receptor]]*Ey
    #   ) %>% as.data.frame

  #** 这里修改了一下，使用qnet进行计算而不是receptor来计算

  n_row <- rep(1,nrow(field_df)) %>% as.matrix() 
  
  field_Fx <- n_row %*% t(ori_field_force$Fx) %>% methods::as("dgCMatrix")
  field_Fy <- n_row %*% t(ori_field_force$Fy) %>% methods::as("dgCMatrix")
  field_Fnorm <- sqrt(field_Fx^2+field_Fy^2) %>% methods::as("dgCMatrix")

  field_force <- 
    list(force_x=field_Fx,force_y=field_Fy,force_norm=field_Fnorm)
  field_force <- 
    lapply(field_force,function(mat){
      rownames(mat) <- colnames(mat) <- rownames(field_df)
      return(mat)
    })
  return(field_force)
}


#' calculation field-spot interaction force
#'
#' calculation field-spot interaction force which 
#' considered as how the whole field affect the spot
#'
#' @param field_df data.frame, contain vector field info.
#' @param ligand Ligand gene.
#' @param receptor Receptor gene.
#' @param LR_mat dgCMatrix, spot2spot interaction matrix, can be LR direcional expression or 
#' LR specific force matrix
#'
#' @return return list of matrices including force component of x,y and force norm.
#' 
#' @importClassesFrom Matrix dgCMatrix CsparseMatrix
#' 
#' @seealso calc_field_force_mat
#' 
#' @export
calc_field_force_SpMat <- function(
  field_df,
  ligand,receptor,
  LR_mat
){
  if(all(rownames(field_df) %in% rownames(LR_mat))){
    if(!all(rownames(field_df) == rownames(LR_mat))){
      field_df <- field_df[rownames(LR_mat),]
    }
  }else(
    stop("rownames of field_df should be same as rownames of LR_mat")
  )

  if(!inherits(LR_mat,"dgCMatrix")){
    LR_mat %<>% methods::as("CsparseMatrix")
  }

  q_net <- field_df[,ligand] - field_df[,receptor]
  Fx=q_net*field_df[,"Ex"]
  Fy=q_net*field_df[,"Ey"]
  
  field_Fx <- field_Fy <- field_Fnorm <- LR_mat

  field_Fx@x <- rep(Fx,times=diff(LR_mat@p)) 
  field_Fy@x <- rep(Fy,times=diff(LR_mat@p)) 
  field_Fnorm@x <- sqrt(field_Fx@x^2+field_Fy@x^2) 

  field_force <- 
    list(force_x=field_Fx,force_y=field_Fy,force_norm=field_Fnorm)
  # field_force <- 
  #   lapply(field_force,function(mat){
  #     rownames(mat) <- colnames(mat) <- rownames(field_df)
  #     return(mat)
  #   })
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
#' @importClassesFrom Matrix dgCMatrix CsparseMatrix
#' @importFrom methods as
#' @export
calc_S2S_score_mat <- function(field_force_list,S2S_force_list){

  field_Fx <- field_force_list$force_x
  field_Fy <- field_force_list$force_y
  field_Fnorm <- field_force_list$force_norm
  
  S2S_Fx <- S2S_force_list$force_x
  S2S_Fy <- S2S_force_list$force_y
  S2S_Fnorm <- S2S_force_list$force_norm

  Ta <- S2S_Fx@x*field_Fx@x
  Tb <- S2S_Fy@x*field_Fy@x
  temp <- (Ta+Tb)


  temp <- temp/field_Fnorm@x
  #S2S_cos <- temp/(field_Fnorm*S2S_Fnorm)
  #df$cos_norm <- (df$Fx*vec$Fx+df$Fy*vec$Fy)/sqrt(vec$Fx^2+vec$Fy^2)

  #** devided by 0 introducing NaN, and change the class into dgeMatrix
  temp[is.na(temp)] <- 0
  S2S_cos_norm <- field_Fnorm
  S2S_cos_norm@x <- temp
  S2S_cos_norm %<>% Matrix::drop0()
  #S2S_cos_norm %<>% methods::as("CsparseMatrix")

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
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
calc_database_S2S_score <- function(kept_db,prep_list,db_S2S_force_list){
  dist_list <- prep_list$dist_list
  LR_mat_list <- prep_list$q_list
  field_df_list <- prep_list$field_list
  
  db_S2S_score_list <- list()
  pb <- utils::txtProgressBar(min = 0, max = nrow(kept_db), style = 3)
  iter=seq_len(nrow(kept_db))
  for(i in iter){
    db_field_force <- 
      calc_field_force_SpMat(
        field_df_list[[i]],
        kept_db$Ligand[i],
        kept_db$Receptor[i],
        db_S2S_force_list[[i]][[1]]
        )
    
    db_S2S_score_list[[i]] <- 
      calc_S2S_score_mat(
        db_field_force,
        db_S2S_force_list[[i]]
        )
    utils::setTxtProgressBar(pb, i)
    
  }
  close(pb)
  names(db_S2S_score_list) <- kept_db$id
  
  return(db_S2S_score_list)
}

#** Cluster2Cluster Interaction Score

#' calculation of cluster2cluster interaction score matrix
#'
#' calculation of cluster2cluster interaction score matrix from sum-up S2S score by clusters
#'
#' @param S2S_score_mat S2S score matrix
#' @param clu_info_list A named list, where each element is a vector. 
#' @param n_mat matrix, each element is the product of number of interacting clusters.
#' Col as sender, row as receiver.
#' The names of the list elements represent the cluster names.
#' Each vector corresponds to a cluster, and the elements within each vector are the indices of the matrix.
#' 
#' @return cluster2cluster interaction score matrix
#' @export
calc_C2C_mat <- function(S2S_score_mat,clu_info_list,n_mat){
  n <- length(clu_info_list)
  C2C_score_mat <- matrix(nrow=n, ncol=n)
  #idx_list <- lapply(clu_info_list, function(clu) which(rownames(S2S_score_mat) %in% clu))
  for (i in 1:n) {
    tmp_idx_i <- clu_info_list[[i]]
    for (j in 1:n) {
      tmp_idx_j <- clu_info_list[[j]]
      C2C_score_mat[i,j] <- 
        sum(S2S_score_mat[tmp_idx_i, tmp_idx_j]) / n_mat[i,j]
    }
  }
  return(C2C_score_mat)
}

#' calculation of cluster2cluster interaction score and p-value
#'
#' calculation of cluster2cluster interaction score and p-value 
#' from sum-up S2S score by clusters and shuffle test for p-value
#'
#' 
#' @inheritParams calc_C2C_mat
#' @inheritParams auto_select_lapply
#' @param clu_shuf_list A named list,similar to clu_info_list, 
#' but each element is a vector of shuffled indices. 
#' 
#' @return list of unfiltered cluster2cluster interaction score summary, score and p-value matrices
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom pbapply pboptions
#' @export
calc_C2C_score_loop <- function(
  S2S_score_mat,clu_info_list,clu_shuf_list,n_mat,
  verbose=TRUE,backend="auto"
){
  #pboptions see '?pbapply::pboptions' example
  # if(verbose){
  #   opb <- pbapply::pboptions(type="none")
  #   on.exit(pbapply::pboptions(opb))
  # }

  #clac the sum force from S2S to C2C
  C2C_score <- calc_C2C_mat(S2S_score_mat,clu_info_list,n_mat)

  C2C_p_value <- list()
  shuffle_iter <- length(clu_shuf_list)
  
  # for(i in 1:shuffle_iter){
  #   C2C_shuf <- calc_C2C_mat(S2S_score_mat,clu_shuf_list[[i]])
  #   #C2C_p_value[C2C_shuf > C2C_score] <- C2C_p_value[C2C_shuf > C2C_score] + 1
  #   C2C_p_value[[i]] <- (C2C_shuf > C2C_score)
  # }

  C2C_p_value <- 
    auto_select_lapply(clu_shuf_list, function(clu_shuf){
      C2C_shuf <- calc_C2C_mat(S2S_score_mat,clu_shuf,n_mat)
      #C2C_p_value[C2C_shuf > C2C_score] <- C2C_p_value[C2C_shuf > C2C_score] + 1
      return(C2C_shuf > C2C_score)
    },verbose = verbose,backend = backend
    )
  C2C_p_value <- Reduce("+", C2C_p_value) / shuffle_iter  

  rownames(C2C_score) <- names(clu_info_list)
  colnames(C2C_score) <- names(clu_info_list)
  rownames(C2C_p_value) <- names(clu_info_list)
  colnames(C2C_p_value) <- names(clu_info_list)

  score_long <- 
  C2C_score %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source,names_to = "Target",values_to = "raw_score") %>%
    dplyr::mutate(scale_score=scale(.data[["raw_score"]])[,1])

  p_value_long <- 
  C2C_p_value %>% 
    as.data.frame() %>%
    tibble::rownames_to_column(var="Source") %>%
    tidyr::pivot_longer(!Source,names_to = "Target",values_to = "p_value")

  long_df <- dplyr::full_join(score_long,p_value_long,by=c("Source","Target"))
  return(list(summary=long_df,score_mat=C2C_score,p_val_mat=C2C_p_value))
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
#' @param shuffle_iter number of shuffle iterations, default is 343, 
#' the cube of 7 and the Spark.
#' @param random_seed random seed used in shuffle. Default is 42, 
#' the answer to the ultimate question of life, the universe, and everything.
#' @inheritParams auto_select_lapply
#' 
#'
#' @return list of each LR pair C2C score summary, score and p value.
#' @export
calc_database_C2C_score <- function(
  kept_db,db_S2S_score_list,
  seurat_obj,cluster=NULL,shuffle_iter=343,random_seed=42,
  verbose=TRUE
){
  cluster_info <- cluster_info_identifier(seurat_obj,cluster)
  #clu_info_list <- split(cluster_info$barcode,cluster_info$cluster)
  clu_info_list <- split(seq_len(nrow(cluster_info)),cluster_info$cluster)
  n_vec <- lapply(clu_info_list,length) %>% as.numeric()
  n_mat <- n_vec %*% t(n_vec)

  set.seed(random_seed)
  clu_shuf_list <- generate_shuffle_list(cluster_info,shuffle_iter)
  #TODO 需要添加一个防呆设置检测输入s2s_force 还是S2S_score

  C2C_score_list <- list()
  iter=seq_len(nrow(kept_db))

  if(length(iter) >= length(clu_shuf_list)){
    backend=c("auto","none")
  }else{
    backend=c("none","auto")
  }

  C2C_score_list <- 
  auto_select_lapply(db_S2S_score_list,function(S2S_score){
      calc_C2C_score_loop(
        S2S_score,
        clu_info_list,
        clu_shuf_list,
        n_mat,
        verbose=FALSE,backend = backend[2]
      )
  },verbose=verbose,backend = backend[1]
  )
  names(C2C_score_list) <- kept_db$id
  set.seed(NULL)

  return(C2C_score_list)
}