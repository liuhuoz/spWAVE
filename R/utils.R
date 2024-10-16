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


#' Concentrate Coordinates and Expression
#' 
#' concentrate expression matrix and coordinates by barcode 
#'
#' @param expr expression matrix.
#' @param genes select genes.
#' @param coord_info spot coordinates, 2 cols with x and y.
#'
#' @return return dataframe include barcode, spot coordinates, and selected genes expression values.

concentrate_coord_expr <- function(expr,genes,coord_info){
  expr_df <- expr[,genes,FALSE]
  check_result <- check_row_order(rownames(expr_df),rownames(coord_info))
  if(check_result){
    merge_df <- cbind(coord_info,expr_df)
  }else{
    expr_df$barcode <- rownames(expr_df)
    merge_df <- left_join(coord_info,expr_df,by="barcode")
    rownames(merge_df) <- merge_df$barcode
    merge_df <- merge_df[rownames(coord_info),]
  }
  return(merge_df)
}



#' Cluster Info Identifier
#'
#' Retrieve cluster info from idents or metadata in Seurat object, or merge the info from given data.frame.
#'
#' @param seurat_obj Seurat object.
#' @param cluster cluster info. Default is NULL and will using active.ident in Seurat object. 
#' Accepted data format is data.frame or character. See detailed for more explain.
#'
#' @return return dataframe include barcode and cluster id.
#' @details The default 'cluster' is NULL and will using active.ident in Seurat object. 
#' If a character is given, it must be the column name in Seurat metadata, and the column will be used as cluster. 
#' If a data.frame is given, the column names should be 'cluster' and 'barcode', respectively.

cluster_info_identifier <- function(seurat_obj,cluster=NULL){
  if(is.null(cluster)){
    cluster_key <- "ident"
  }else if(is.data.frame(cluster)){
    cluster_key <- "df"
  }else if(is.character(cluster)){
    cluster_key <- "char"
  }else{
    stop("Invalid cluster argument")
  }

  switch(cluster_key,
    ident={
      cluster_df <- seurat_obj@active.ident
      cluster_df %<>% as.data.frame()
      cluster_df$barcode <- rownames(cluster_df)
      colnames(cluster_df) <- c("cluster","barcode")
    },

    df={
      res1 <- "cluster" %in% colnames(cluster)
      res2 <- "barcode" %in% colnames(cluster)
      if(res1 & res2){
          cluster_df <- cluster
        }else{
          stop("Invalid data.frame format, check colnames")
        }
    },

    char={
      if(cluster %in% colnames(seurat_obj@meta.data)){
        meta <- seurat_obj@meta.data
        meta$barcode <- rownames(meta)
        cluster_df <- meta[,c("barcode",cluster)]
        colnames(cluster_df) <- c("barcode","cluster")
      }else{
        stop("Cluster not detected")
      }
    }
  )

  return(cluster_df)
}


#' Concentrate LR expression and field estimate result
#'
#' Concentrate LR expression and field estimate result for each spot or cell
#'
#' @param db_field_result field estimate result of LR in database .
#' @param ligand Ligand genes.
#' @param receptor Receptor genes or complex.
#'
#' @return return dataframe include barcode, spot coordinates, genes and complex expression and field estimate result.

concentrate_LR_field_info <- function(db_field_result,ligand,receptor){
  single_mol_field_list <- db_field_result$single_mol_field_list
  LR_pair_field_list <- db_field_result$LR_pair_field_list
  LR_pair <- paste(ligand,receptor,sep=".")

  merge_df <- cbind(
      LR_pair_field_list[[LR_pair]],
      single_mol_field_list[[ligand]][,ligand,FALSE],
      single_mol_field_list[[receptor]][,receptor,FALSE]
      )
  return(merge_df)
}


#' Generate Shuffle Cluster List
#'
#' Generate shuffle list of barcode and cluster,for downstream permutation test
#'
#' @param cluster_info a dataframe with columns "barcode" and "cluster"
#' @param shuffle_iter Integer, number of shuffle iterations
#'
#' @return List of dataframe, each dataframe contains shuffled barcode and cluster
#' and the length of list will be equal to shuffle_iter

generate_shuffle_list <- function(cluster_info,shuffle_iter=500){
  clu_shuf_list <- list()
  for (i in 1:shuffle_iter){
    clu_shuffle <- 
    data.frame(
      index=sample(1:nrow(cluster_info)),
      cluster=cluster_info$cluster
    )
    clu_shuf_list[[i]] <- split(clu_shuffle$index,clu_shuffle$cluster)
  }
  names(clu_shuf_list) <- seq_len(shuffle_iter)
  return(clu_shuf_list)
}


#' Aggregate C2C_score form list 
#'
#' Aggregate C2C_score summaries from field estimate results list, for plotting.
#'
#' @param db_C2C_score_list The list of C2C_score field estimate results. 
#' @param kept_db filtered database, should be match the C2C_score_list.
#' 
#' @details This function extract C2C_score summary dataframe from the list of C2C_score field estimate results.
#' And do not filter anything and all value will be kept, including the no significant p value and minus value of C2C_score which often considered as reverse signal direcetion.
#' 
#' @return return a dataframe contain cluster id and C2C_score of each LR pair or family.
aggregate_C2C_score <- function(db_C2C_score_list,kept_db){
  db_prep_list <- lapply(db_C2C_score_list,function(x) x$summary)
  db_prep_list <- 
    lapply(db_prep_list,function(db){
      db %<>% 
        #filter(p_value<0.05) %>%
        arrange(desc(raw_score))
      return(db)
    })

  iter <- seq_len(length(db_C2C_score_list))
  for(i in iter){
    db_prep_list[[i]]$id <- names(db_C2C_score_list)[i]
  }
  aggregate_df <- do.call(rbind,db_prep_list)

  aggregate_df %<>% 
    arrange(desc(raw_score)) %>% 
    left_join(kept_db[,c("id","Family")],by="id")
  return(aggregate_df)
} #** for plot



#' Prepare C2C score dataframe for plotting
#'
#' Filter C2C score dataframe by certain keyword and pre-process the data for plotting
#'
#' @param C2C_score_df aggregated C2C score dataframe from aggregate_C2C_score. Contain cluster id and score of each LR pair or family.
#' @param LR_pair filter LR pair.
#' @param LR_family filter LR family.
#' @param source_use filter cluster id as source (sender).
#' @param target_use filter cluster id as target (receiver).
#' @param scale scale the C2C score. Default is TRUE.
#'
#' @return filtered C2C score dataframe.

prep_C2C_plot_df <- function(
  C2C_score_df,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  scale=TRUE
){
  LR_pair <- filter_cat_keyword(C2C_score_df,"id",LR_pair)
  LR_family <- filter_cat_keyword(C2C_score_df,"Family",LR_family)
  source_use <- filter_cat_keyword(C2C_score_df,"Source",source_use)
  target_use <- filter_cat_keyword(C2C_score_df,"Target",target_use)

  filtered_df <- 
    C2C_score_df %>%
      filter(
        id %in% LR_pair,
        Family %in% LR_family,
        Source %in% source_use,
        Target %in% target_use
      )
  
  if(scale){
    filtered_df$scale_score <- 
      scale(filtered_df$raw_score,center=TRUE,scale=TRUE)
  }

  return(filtered_df)
}

#' Filter category by keyword
#' 
#' Filter dataframe by category and keyword
#' 
#' @param df Dataframe to be filtered.
#' @param filter_cat The category to be filtered. Must be a colname of the dataframe.
#' @param filter_key the keyword to be filtered, should be found in the category. 
#' Default is NULL. If NULL, all result will be kept.
#'
#' @return Dataframe filtered by category and respective keyword.

filter_cat_keyword <- function(df,filter_cat,filter_key=NULL){
  all_key <- table(df[[filter_cat]]) %>% names
  if(is.null(filter_key)){
    filter_key <- all_key
  }else if(!all(filter_key %in% all_key)){
    stop(paste0("Filter keywords not found in ",filter_cat))
  }
  return(filter_key)
}



#' Assign colors to cluster for chord plot
#'
#' Assign color to cluster for chord plot with built-in Material Design color palette.
#' 
#' @param C2C_score_df C2C score dataframe contain cluster in column "Source".
#' 
#' @details The color will be assigned by the order of "Source" and "Target" seprately in C2C_score_df, and specifically used in the chord plot. 
#' Typically, the input dataframe is unfiltered and contains all cluster in "Source". For more information, please refer to \code{\link{aggregate_C2C_score}}
#' To keep the color consistent crossing the clusters, this function only using Source to assign the color.
#' This function return 2 side result of cluster labelled by Source and Target.
#' For more universal color assignment, please use \code{\link{assign_clu_col_lite}}.
#' 
#' @return return a named and ordered vector of colors of both Source and Target, seprately. The names are cluster id with prefix "S@" and "R@".
#' 
#' @examples
#'#' C2C_score_df <- 
#     aggregate_C2C_score(db_C2C_score_list,kept_db)
#' grid_col <- assign_cluster_color(C2C_score_df)
#' 
assign_cluster_color <- function(C2C_score_df){
  pic_df <- 
    C2C_score_df %>%
      #filter(p_value<0.05) %>%
      #filter(Source!=Target) %>%
      mutate(Source = paste0("S@",Source),Target = paste0("R@",Target))

  clu_name <- 
    C2C_score_df$Source %>%
    table() %>% names() 

  chord_order <- 
    c(
      paste0(rep("S@",length(clu_name)),clu_name),
      paste0(rep("R@",length(clu_name)),rev(clu_name))
    )
  grid_col <- 
    c(MD2_color_picker(length(clu_name)), 
      rev(MD2_color_picker(length(clu_name))))
  names(grid_col) <- chord_order

  temp <- which(chord_order %in% c(pic_df$Source,pic_df$Target))
  chord_order <- chord_order[temp]
  grid_col <- grid_col[temp]

  return(grid_col)
}


#' Universal assign colors to cluster for plot
#'
#' Assign color to cluster for plot with built-in Material Design color palette.
#'
#' @param C2C_score_df C2C score dataframe contain cluster in column "Source".
#' 
#' @details The color will be assigned by the order of "Source" in C2C_score_df.
#' Typically, the input dataframe is unfiltered and contains all cluster in "Source". For more information, please refer to \code{\link{aggregate_C2C_score}}
#' To keep the color consistent crossing the clusters, this function only using Source to assign the color
#' This function only return one side result of cluster which will not distinguish the Source or Target.
#' 
#' @return return a named vector of colors, the names are cluster id.
#' 
#' @examples
#'#' C2C_score_df <- 
#     aggregate_C2C_score(db_C2C_score_list,kept_db)
#' grid_col <- assign_clu_col_lite(C2C_score_df)
assign_clu_col_lite <- function(C2C_score_df){
  clu_name <- 
    C2C_score_df$Source %>%
    table() %>% names() 

  clu_col <- MD2_color_picker(length(clu_name))
  names(clu_col) <- clu_name

  return(clu_col)
}