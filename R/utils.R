#' get_spatial_expr
#'
#' Exract coordinates and expression values from Seurat object
#'
#' @param SrtObj Seurat object.
#' @param gene one or more gene names in Seurat object.
#' @param assay assay name in slot of Seurat object. Default="SCT"
#'
#' @return return a data.frame containing coordinates, barocdes, and expression values.
#' 
#' @details Right now this function only test on Xenium and Visium. 
#' The unit of return coordinate is micrometer after scaling. 
#' For scale factor of Visium, please see \url{https://github.com/satijalab/seurat/issues/6980} \cr
#' To keep the consistent between plot systems of Seurat and directly using ggplot2,
#' the coordiantes were rotated differently based on the ST techniques.
#' For more infos, please see: \cr
#' Visium rotation: \url{https://github.com/satijalab/seurat/issues/2702} \cr
#' Xenium rotation: \url{https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788}
#' @importFrom Seurat GetTissueCoordinates GetAssayData
#' @import magrittr
#' @import dplyr
#' @export
get_spatial_expr <- function(SrtObj,gene,assay="SCT"){
  if(class(SrtObj@images[[1]]) %in% c("VisiumV1","VisiumV2")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj,scale=NULL)
    scale_spot <- SrtObj@images[[1]]@scale.factors$spot
    scale_micrometer <- 65/scale_spot
    coord_info$imagerow <-
      max(coord_info$imagerow) - coord_info$imagerow + min(coord_info$imagerow)
    coord_info <- coord_info * scale_micrometer
    coord_info$barcode <- rownames(coord_info)
    colnames(coord_info) <- c("y","x","barcode")
    #** scale the image coordinate unit pixel to real world length unit micrometer
    #** for more info, see https://github.com/satijalab/seurat/issues/6980
    #** rotation of coordinate because of the plot consistent using ggplot2 and Seurat plot
    #** for more info, see https://github.com/satijalab/seurat/issues/2702
  }else if(class(SrtObj@images[[1]]) %in% c("FOV")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj)
    colnames(coord_info) <- c("y","x","barcode")
    #** here we exchange x and y as coord_flip, which was introduced by Seurat
    #** which will refine the plot direction in ggplot2.
    #** for more info, see https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788
  }

  expr_mat <- Seurat::GetAssayData(SrtObj,assay=assay)
  select_gene <-
    expr_mat %>% as.matrix() %>% t() %>%
    as.data.frame() %>%
    dplyr::select(all_of(gene))
  colnames(select_gene) <- gene
  select_gene$barcode <- rownames(select_gene)

  spatial_expr <- dplyr::full_join(coord_info,select_gene,by="barcode")
  #colnames(spatial_expr) <- c("y","x","barcode",gene)
  return(spatial_expr)
}

#' Get coordinates from Seurat Object
#'
#' Extract coordinates from Seurat object, and modify for calculation
#'
#' @param SrtObj Seurat object.
#'
#' @return return a data.frame containing coordinates and barcode or cell id.
#' 
#' @details Right now this function only test on Xenium and Visium. 
#' The unit of return coordinate is micrometer after scaling. 
#' For scale factor of Visium, please see \url{https://github.com/satijalab/seurat/issues/6980} \cr
#' To keep the consistent between plot systems of Seurat and directly using ggplot2,
#' the coordiantes were rotated differently based on the ST techniques.
#' For more infos, please see: \cr
#' Visium rotation: \url{https://github.com/satijalab/seurat/issues/2702} \cr
#' Xenium rotation: \url{https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788}
#' @importFrom Seurat GetTissueCoordinates
#' @export
get_coordinates <- function(SrtObj){
  if(class(SrtObj@images[[1]]) %in% c("VisiumV1","VisiumV2")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj,scale=NULL)
    scale_spot <- SrtObj@images[[1]]@scale.factors$spot
    scale_micrometer <- 65/scale_spot
    coord_info$imagerow <-
      max(coord_info$imagerow) - coord_info$imagerow + min(coord_info$imagerow)
    coord_info <- coord_info * scale_micrometer
    coord_info$barcode <- rownames(coord_info)
    colnames(coord_info) <- c("y","x","barcode")
    #** scale the image coordinate unit pixel to real world length unit micrometer
    #** for more info, see https://github.com/satijalab/seurat/issues/6980
    #** rotation of coordinate because of the plot consistent using ggplot2 and Seurat plot
    #** for more info, see https://github.com/satijalab/seurat/issues/2702
  }else if(class(SrtObj@images[[1]]) %in% c("FOV")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj)
    colnames(coord_info) <- c("y","x","barcode")
    rownames(coord_info) <- coord_info$barcode
    #** here we exchange x and y as coord_flip, which was introduced by Seurat
    #** which will refine the plot direction in ggplot2.
    #** for more info, see https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788
  }

  return(coord_info)
}

#' spa_vectorized_pdist
#' 
#' @details from an excellent post: https://www.r-bloggers.com/2013/05/pairwise-distances-in-r/
#' this function is called by other functions to quickly compute the distance between
#' cells to grid points, or between grid points
#'
#' @param source matrix
#' @param target matrix
#' @return returns pairwise-distances
#' 
#' @export
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
#' @details pairwise calculate 2 matrix subtraction by row between
#' cells to grid points, or between grid points,and return a list of each target row.
#' The result would be source subtracted by target.
#'
#' @param source matrix
#' @param target matrix
#' @return returns pairwise-subtract for each row list
#'  
#' @export
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
#' cpdb_human <- load_database(db_source="CellPhoneDB",db_species="human")
#' @export 
load_database <- function(
  db_source=c("CellChat","CellPhoneDB"),
  db_species=c("human","mouse","zebrafish"),
  filter_type="None"
){
  if(db_source=="CellPhoneDB" && db_species=="zebrafish"){
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
#' @importFrom Seurat GetAssayData
#' @import stringr
#' @return filtered expression matrix
#' @export 
filter_LR_expr <- function(
  db,
  SrtObj,
  assay="RNA",
  min_expr=0.1,
  min_n_cell=NULL,
  min_pct_cell=0.01
){
  expr_mat <- Seurat::GetAssayData(SrtObj,assay=assay)
  db_all_genes <- 
    c(db$Ligand,stringr::str_split(db$Receptor,pattern = "_")) %>% 
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
#' @export 
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
#' @param db dataframe of LR database.
#' @param expr_mat gene expression matrix.
#' @param complex_min_cell minimum number of cell expressing complex
#' 
#' @details The receptor complex expression are calculated using 
#' \deqn{R_{complex} = \prod_{i \in R} R_i ^{1/n_{R}}} where 
#' R is the expression of each sub-unit and n_{R} is the number of sub-units in
#' the complex.
#' 
#' @return Ruturn a list including filtered database, LR expression list, and merged LR expression dataframe.
#' @importFrom matrixStats rowProds
#' @import stringr
#' @export 
generate_complex_data <- function(db,expr_mat,complex_min_cell=10){
  expr_df <- expr_mat %>% as.matrix() %>% t() %>% as.data.frame()
  expr_LR_list <- list()
  kept_db_index <- c()
  for(i in seq_len(nrow(db))){
    lig <- db$Ligand[i]
    rec_complex <- db$Receptor[i]
    rec <- stringr::str_split_1(rec_complex,pattern = "_")
    if(all(c(lig,rec) %in% colnames(expr_df))){ #filter out not-exist in expr LR_pair
      LR_complex_df <- expr_df[,c(lig,rec)]
      temp <- matrixStats::rowProds(LR_complex_df[,c(rec)] %>% as.matrix())^(1/length(rec))
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
  return(list(kept_db=kept_db,
              #expr_LR_list=expr_LR_list,
              expr_LR_df=expr_LR_df))
}
#** 发现expr_LR_list中会存在colnames 重复的情况，不过考虑到一个部分后面不会用到，先不考虑特别处理。
#** 甚至后续可以把这一项删除。
#' Filter out subunit 
#'
#' Filter out subunit 
#'
#' @param kept_db dataframe of LR database.
#' @param expr expression matrix.
#' 
#' @return expr dataframe without subunit expression
filter_subunit <- function(kept_db,expr){
    all_LR <- c(kept_db$Ligand,kept_db$Receptor) %>% unique()
    kept_LR <- intersect(all_LR,colnames(expr))
    return(expr[,kept_LR])
}

#' concatenate Coordinates and Expression
#' 
#' concatenate expression matrix and coordinates by barcode 
#'
#' @param expr expression matrix.
#' @param genes select genes.
#' @param coord_info spot coordinates, 2 cols with x and y.
#'
#' @return return dataframe include barcode, spot coordinates, and selected genes expression values.
#' @export 
concatenate_coord_expr <- function(expr,genes,coord_info){
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


#' Check Rownames Order Consistent
#'
#' Check rownaems are same with standard rownames
#'
#' @param check_row rownames to be checked.
#' @param standrad_row rownames as standard.
#'
#' @return logical TRUE or FALSE
#' @export
check_row_order <- function(check_row,standrad_row){
  row_result <- ifelse(check_row==standrad_row,TRUE,FALSE)
  check_result <- all(row_result)
  return(check_result)
}

#' Cluster Info Identifier
#'
#' Retrieve cluster info from idents or metadata in Seurat object, or merge the info from given data.frame.
#'
#' @param seurat_obj Seurat object.
#' @param cluster cluster info. Default is NULL and will using active.ident in Seurat object. 
#' Accepted data format is data.frame or character. See details for more explain.
#'
#' @return return dataframe include barcode and cluster id.
#' @details The default 'cluster' is NULL and will using active.ident in Seurat object. 
#' If a character is given, it must be the column name in Seurat metadata, and the column will be used as cluster. 
#' If a data.frame is given, the column names should be 'cluster' and 'barcode', respectively.
#' @export 
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


#' concatenate LR expression and field estimate result
#'
#' concatenate LR expression and field estimate result for each spot or cell
#'
#' @param db_field_result field estimate result of LR in database .
#' @param ligand Ligand genes.
#' @param receptor Receptor genes or complex.
#'
#' @return return dataframe include barcode, spot coordinates, genes and complex expression and field estimate result.
#' @export 
concatenate_LR_field_info <- function(db_field_result,ligand,receptor){
  single_mol_field_list <- db_field_result$single_mol_field_list
  LR_pair_field_list <- db_field_result$LR_pair_field_list
  LR_pair <- paste(ligand,receptor,sep=".")

  merge_df <- cbind(
      LR_pair_field_list[,,LR_pair],
      single_mol_field_list[,"expression",ligand],
      single_mol_field_list[,"expression",receptor]
      )
  colnames(merge_df)[4:5] <- c(ligand, receptor)
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
#' @export 
generate_shuffle_list <- function(cluster_info,shuffle_iter=500){
  clu_shuf_list <- list()
  for (i in 1:shuffle_iter){
    clu_shuffle <- 
    data.frame(
      index=sample(seq_len(nrow(cluster_info))),
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
#' @import dplyr
#' @return return a dataframe contain cluster id and C2C_score of each LR pair or family.
#' @export 

aggregate_C2C_score <- function(db_C2C_score_list,kept_db){
  db_prep_list <- lapply(db_C2C_score_list,function(x) x$summary)
  db_prep_list <- 
    lapply(db_prep_list,function(db){
      db %<>% 
        #dplyr::filter(p_value<0.05) %>%
        dplyr::arrange(dplyr::desc(raw_score))
      return(db)
    })

  iter <- seq_len(length(db_C2C_score_list))
  for(i in iter){
    db_prep_list[[i]]$id <- names(db_C2C_score_list)[i]
  }
  aggregate_df <- do.call(rbind,db_prep_list)

  aggregate_df %<>% 
    dplyr::arrange(desc(raw_score)) %>% 
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
#' @import dplyr
#'
#' @return filtered C2C score dataframe.
#' @export 
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
      dplyr::filter(
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
#' @export 
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
#' Typically, the input dataframe is unfiltered and contains all cluster in "Source". For more information, please refer to \code{\link{aggregate_C2C_score}}.
#' To keep the color consistent crossing the clusters, this function only using Source to assign the color.
#' This function return 2 side result of cluster labelled by Source and Target.
#' For more universal color assignment, please use \code{\link{assign_clu_col_lite}}.
#' 
#' @import dplyr
#' @return return a named and ordered vector of colors of both Source and Target, seprately. The names are cluster id with prefix "S@" and "R@". 
#' @export
assign_cluster_color <- function(C2C_score_df){
  pic_df <- 
    C2C_score_df %>%
      #dplyr::filter(p_value<0.05) %>%
      #dplyr::filter(Source!=Target) %>%
      dplyr::mutate(Source = paste0("S@",Source),Target = paste0("R@",Target))

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
#' To keep the color consistent crossing the clusters, this function only using Source to assign the color.
#' This function only return one side result of cluster which will not distinguish the Source or Target.
#' 
#' @return return a named vector of colors, the names are cluster id.
#' 
#' @export 
assign_clu_col_lite <- function(C2C_score_df){
  clu_name <- 
    C2C_score_df$Source %>%
    table() %>% names() 

  clu_col <- MD2_color_picker(length(clu_name))
  names(clu_col) <- clu_name

  return(clu_col)
}


#' Check color matching cluster info
#'
#' Check the given color is matching with cluster info
#'
#' @param color named vector, the values are colors and the names are cluster.
#' @param C2C_score_df dataframe contain cluster, score, and LR pair.
#' @param source_use the cluster will be used as source, which should be kept in return
#' @param target_use the cluster will be used as target, which should be kept in return
#' 
#' @return return a named filtered color vector.
#' @export 
check_cluster_color <- function(color,C2C_score_df,source_use,target_use){
  all_cluster <- c(C2C_score_df$Source,C2C_score_df$Target) %>% unique()
  if(is.null(names(color))){
    stop("cluster color should be a named vector")
  }else if(!all(names(color) %in% c(source_use,target_use,all_cluster))){
    stop("cluster color names should be same as cluster names")
  }
  if(any(duplicated(color))){
    warning("colors are not unique")
  }
  return(color)
}


#' Prepare matrix for heatmap and dot plot
#'
#' Filter and pre-process matrix for heatmap and dot plot
#'
#' @param C2C_score_df C2C score dataframe contain cluster, score, and LR pair.
#' @param LR_pair filter LR pair display in plot.
#' @param LR_family filter LR family display in plot.
#' @param source_use filter cluster as `Source`
#' @param target_use filter cluster as `Target`
#' @param method_use Statistical method to aggregate the C2C score.
#' Only support "count" and "strength". Default is "count".
#' The 'count' reflect the number of interactions,
#' and 'strength' reflect the sum of C2C score.
#' 
#' @import dplyr
#' @importFrom tidyr pivot_wider
#'
#' @return return a filtered matrix which column as Target and row as Source,
#' and the value based on `method_use`.
#' @export 
prep_plot_matrix <- function(
  C2C_score_df,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  method_use=c("count","strength")
){
  ht_df <- 
    C2C_score_df %>%
      prep_C2C_plot_df( 
        LR_pair=LR_pair,
        LR_family=LR_family,
        source_use=source_use,
        target_use=target_use,
        scale=FALSE
      ) %>%
      dplyr::filter(p_value<0.05)

  method_use=match.arg(method_use)
  method_use <-ifelse(method_use=="count","length","sum")

  ht_mat <- 
  ht_df[,c(1:3)] %>%
    tidyr::pivot_wider(
      id_cols = Source,
      names_from = Target,
      values_from = raw_score,
      values_fn = list(raw_score = method_use)) %>%
      as.data.frame()
  ht_mat[is.na(ht_mat)] <- 0

  rownames(ht_mat) <- ht_mat$Source
  ht_mat <- ht_mat[,-1,FALSE] %>% as.matrix()

  #** extend matrix to include non-filtered but 0 interaction cluster
  kept_row <- filter_cat_keyword(C2C_score_df,"Source",source_use)
  kept_col <- filter_cat_keyword(C2C_score_df,"Target",target_use)

  temp <- which(!(kept_row %in% rownames(ht_mat)))
  temp_mat <- matrix(0,nrow = length(temp),ncol = ncol(ht_mat))
  rownames(temp_mat) <- kept_row[temp]
  colnames(temp_mat) <- colnames(ht_mat)
  ht_mat <- rbind(ht_mat,temp_mat)

  temp <- which(!(kept_col %in% colnames(ht_mat)))
  temp_mat <- matrix(0,nrow = nrow(ht_mat),ncol = length(temp))
  rownames(temp_mat) <- rownames(ht_mat)
  colnames(temp_mat) <- kept_col[temp]
  ht_mat <- cbind(ht_mat,temp_mat)

  ht_mat <- ht_mat[kept_row,kept_col]
  return(ht_mat)
}



#' Extract LR field result form results list
#'
#' Extract certain LR field form complex list and aggregate into a dataframe
#'
#' @param database_result database result list from spWAVE.
#' @param LR LR pair or family. LR pair must be format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param kept_db The database used in results list.
#' 
#' @import dplyr
#' @import abind
#'
#' @return return a dataframe contain cluster, coordinates, 
#' expression of each ligand and receptor and field estimate result including direction and strength.
#' @export 
extract_LR_field_result <- function(database_result,LR,kept_db){
  db_res <- database_result
  LR_pool <- abind::abind(
    db_res$LR_pair_field_list, 
    db_res$LR_family_field_list, along = 3)

  if(!(LR %in% dimnames(LR_pool)[[3]])){
    stop("LR not found!")
  }else{
    LR_field <- LR_pool[,,LR]
  }

  #** 分开获取LR的名字，方便后面进行加减
  L_info <- 
    kept_db %>% 
      dplyr::filter(Family==LR | id ==LR) %>%
      dplyr::pull(Ligand) %>%
      unique()
  R_info <- 
    kept_db %>% 
      dplyr::filter(Family==LR | id ==LR) %>%
      dplyr::pull(Receptor) %>%
      unique()

  #** 数据库计算的单分子向量场结果df内，基因表达量都在第4列，所以这里直接提取第4列就行
  #** 使用array后再次修改,直接进行切片操作，丢弃维度2即可。
  LR_expr_df <- 
    db_res$single_mol_field_list[,"expression",c(L_info,R_info),drop=TRUE]

  df <- cbind.data.frame(LR_field,LR_expr_df)
  df %<>% 
    mutate(Rel_LR_Exp=
      rowSums(across(all_of(L_info)))-
      rowSums(across(all_of(R_info)))
    )
  return(df)
}



#' optimize arrow size for field projection
#'
#' optimize arrow size for plot_field_direction function
#'
#' @param grid.info coordinates and vector (Ex,Ey) for each point
#' @param scale.factor sacle.factor, smaller value means larger size.
#' @param normalize logical default is FALSE. True means all vector length will be 1
#'
#' @return return a data.frame with optimize arrow size for plotting.
#' @export
optimize.arrow <- function(grid.info, scale.factor = 1,normalize=FALSE){
  grid.info$qsum <- sqrt(grid.info$Ex^2 + grid.info$Ey^2)
  dim.scale <- mean(max(grid.info$x) - min(grid.info$x), 
                    max(grid.info$y) - min(grid.info$y))
  #for normalized arrow length using this one, and proper scale.factor=50
  #nor <- match.arg(normalize,choices=c(FALSE,TRUE))
  if(normalize){
      grid.info$sf = (dim.scale/scale.factor)/grid.info$qsum
    }else{
      grid.info$sf = (dim.scale/scale.factor)
    }
  grid.info$Ex.u <- grid.info$Ex * grid.info$sf
  grid.info$Ey.u <- grid.info$Ey * grid.info$sf
  return(grid.info)
}


#' Auto calculate arrow scale factor
#'
#' Auto calculate arrow scale factor for plot_field_direction2 function,
#' based on the average magnitude of vector.
#'
#' @param grid_info coordinates and vector (Ex,Ey) for each point
#' @details This function is auto calculated scale factor for plot_field_direction2
#' by average magnitude of vector. 
#' The calculated scale factor is used in "point" mode to control arrow size.
#'
#' @importFrom stats median
#' @return return a number indicating proper scale factor for plotting.
#' @export
autocalc_arrow_sf <- function(grid_info){
  temp <- c(grid_info$Ex,grid_info$Ey) 
  num_power <- temp %>% abs() %>% log10() %>% stats::median() %>% round()
  sf <- 10^(1 - num_power)
  return(sf)
}

#' Get slide aspect ratio
#'
#' calculate the aspect ratio of a subset image from Seurat object
#'
#' @param SeuObj subset Seurat object.
#' 
#' @importFrom Seurat GetTissueCoordinates GetAssayData
#' @return return aspect ratio 
#' 
#' @export 
subset_ratio <- function(SeuObj){
  coord <- Seurat::GetTissueCoordinates(SeuObj)
  # calculate the aspect ratio of rows to columns
  ratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  # force the image into the right aspect ratio
  return(ratio)
}



#' Get coordinates from Seurat Object for plot
#'
#' Extract coordinates from Seurat object, and modify for plotting
#'
#' @param SrtObj Seurat object.
#'
#' @return return a data.frame containing coordinates and barcode or cell id.
#' 
#' @details Right now this function only test on Xenium and Visium. 
#' The unit of return coordinate is micrometer after scaling. 
#' For scale factor of Visium, please see \url{https://github.com/satijalab/seurat/issues/6980} \cr
#' To keep the consistent between plot systems of Seurat and directly using ggplot2,
#' the coordiantes were rotated differently based on the ST techniques.
#' For more infos, please see: \cr
#' Visium rotation: \url{https://github.com/satijalab/seurat/issues/2702} \cr
#' Xenium rotation: \url{https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788}
#' @importFrom Seurat GetTissueCoordinates
#' @export
get_coordinates_in_plot <- function(SrtObj){
  if(class(SrtObj@images[[1]]) %in% c("VisiumV1","VisiumV2")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj)
    scale_spot <- SrtObj@images[[1]]@scale.factors$spot
    #scale_micrometer <- 65/scale_spot
    coord_info$imagerow <-
      max(coord_info$imagerow) - coord_info$imagerow + min(coord_info$imagerow)
    #coord_info <- coord_info * scale_micrometer
    coord_info$barcode <- rownames(coord_info)
    colnames(coord_info) <- c("y","x","barcode")
    #** scale the image coordinate unit pixel to real world length unit micrometer
    #** for more info, see https://github.com/satijalab/seurat/issues/6980
    #** rotation of coordinate because of the plot consistent using ggplot2 and Seurat plot
    #** for more info, see https://github.com/satijalab/seurat/issues/2702
  }else if(class(SrtObj@images[[1]]) %in% c("FOV")){
    coord_info <- Seurat::GetTissueCoordinates(SrtObj)
    colnames(coord_info) <- c("y","x","barcode")
    rownames(coord_info) <- coord_info$barcode
    #** here we exchange x and y as coord_flip, which was introduced by Seurat
    #** which will refine the plot direction in ggplot2.
    #** for more info, see https://github.com/satijalab/seurat/issues/6110#issuecomment-1172650788
  }

  return(coord_info)
}

#' check character is color
#'
#' check character is color
#'
#' @param x character
#' @details adapted from \url{https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation}
#' @return logical
#' @export 
is.color <- function(x){
  sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
              error = function(e) FALSE)
      })
}



#' generate vector field in grid
#'
#' calculation vector field in grid from vector field in spatial coordinate 
#'
#' @param spatial_vector data.frame,spatial vector field in point coordinate,
#' Must contain x,y,Ex,Ey.
#' @param grid_density grid density, default is 1. Smaller value means sparser grid.
#' @param grid_knn knn number, default is NULL. If NULL, it will be set to 1/50 of the number of points.
#' @param grid_scale grid scale factor, default is 1. Relate to filtering grid point. 
#' @param grid_thresh grid threshold, default is 1. Relate to filtering grid point.
#' @details This part adapted from COMMOT plot_cell_signaling function in R
#' See \url{https://github.com/zcang/COMMOT}
#' @importFrom FNN get.knnx
#' @importFrom stats dnorm quantile
#' @return data.frame, contain x,y,Ex,Ey.
#' @export

generate_grid_vector <- function(spatial_vector, grid_density = 1,grid_knn=NULL,grid_scale = 1.0,grid_thresh = 1.0){
  #This part adapted from COMMOT plot cell_signaling function in R
  #generate squre grid
  X <- spatial_vector[,c("x","y")]
  V <- spatial_vector[,c("Ex","Ey")]

  xl <- min(X[, "x"])
  xr <- max(X[, "x"])
  epsilon <- 0.02 * (xr - xl)
  xl <- xl - epsilon
  xr <- xr + epsilon
  yl <- min(X[, "y"])
  yr <- max(X[, "y"])
  epsilon <- 0.02 * (yr - yl)
  yl <- yl - epsilon
  yr <- yr + epsilon
  ngrid_x <- as.integer(50 * grid_density)
  gridsize <- (xr - xl) / ngrid_x
  ngrid_y <- as.integer((yr - yl) / gridsize)
  meshgrid <- expand.grid(x = seq(xl, xr, length.out = ngrid_x), y = seq(yl, yr, length.out = ngrid_y))
  grid_pts <- cbind(meshgrid$x, meshgrid$y)
  colnames(grid_pts) <- c("x","y")

  #knn filter out outliers
  if (is.null(grid_knn)) {
      grid_knn <- as.integer(nrow(X) / 50)
  }
  nn_mdl <- FNN::get.knnx(X[,1:2], grid_pts,algorithm = "kd_tree",k = grid_knn)
  dis <- nn_mdl$nn.dist
  nbs <- nn_mdl$nn.index
  w <- stats::dnorm(x = dis, mean = 0, sd = gridsize * grid_scale)
  w_sum <- rowSums(w)

  grid_thresh <- grid_thresh * quantile(w_sum, probs = 0.99) / 100
  grid_pts <- grid_pts[w_sum > grid_thresh, ]

  V_grid <- data.frame()
  for(i in seq_len(nrow(nbs))){
    temp <- (V[nbs[i,],]*w[i,]) %>% colSums()
    V_grid <- rbind.data.frame(V_grid,temp)
  }
  colnames(V_grid) <- c("Ex","Ey")
  V_grid <- V_grid / pmax(1, w_sum)
  V_grid <- V_grid[w_sum > grid_thresh,]

  merge_df <- cbind(grid_pts,V_grid)
  return(merge_df)
}


#** for xenium like centroid method

#' The implation of Generate Meta coordinates 
#'
#' The implation of Generate meta coordinates by kmeans methods
#'
#' @param coord data.frame,spatial coordinates, 
#' should have column names "x" and "y" and "barcode".
#' @param random_seed random seed of kmeans for reproducibility. Default is 42.
#' @inheritParams stats::kmeans
#' 
#' @seealso \code{\link[stats]{kmeans}}
#' 
#' @import dplyr
#' @importFrom stats kmeans
#'
#' @return list, meta coordinates contain geometric center of each cluster and 
#' spot ids. And km_cluster contain the barcode and cluster id.
#' @export 
generate_kmeans_coord_impl <- function(
  coord,
  centers=nrow(coord)/10,
  iter.max=10,
  nstart=1,
  random_seed=42
){
  set.seed(random_seed)
  km_result <- stats::kmeans(coord[,c("x","y")], centers=centers, nstart = nstart,iter.max = iter.max)
  km_coord <- 
    cbind.data.frame(
      km_result$centers,
      barcode=paste0("clu_",rownames(km_result$centers))
    )
  coord$cluster <- paste0("clu_",km_result$cluster)
  cluster_info <- coord[,c("barcode","cluster")]
  rownames(km_coord) <- km_coord$barcode
  
  set.seed(NULL)

  return(list(meta=km_coord,km_cluster=cluster_info))
}

#' Generate ROI barcodes
#' 
#' Extract ROI barcodes from Seurat FOV object and extend adjacent barcodes in same clusters.
#' 
#' @param seurat_fov Seurat FOV object, e.g. seurat_obj\[\["fov"\]\].
#' @param clu_info data.frame, cluster info, the colnames must be "barcode" and "cluster".
#' Usally generated by \code{\link{generate_kmeans_coord}}
#' 
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom methods is
#' 
#' @export
generate_ROI_barcode <- function(seurat_fov,clu_info){
  if(!is(seurat_fov,"FOV")){
    stop("Invalid seurat fov object, please check the input.")
  }
  ROI_barcode <- GetTissueCoordinates(seurat_fov)[,3]

  temp <- which(clu_info$barcode %in% ROI_barcode)
  inner_clu <- clu_info$cluster[temp] %>% unique()
  inner_clu_idx <- which(clu_info$cluster %in% inner_clu)
  updated_barcode <- clu_info$barcode[inner_clu_idx]
  return(updated_barcode)
}

#' Generate Holed coordinates
#'
#' Generate a coordinates with single coordinates in a circle,
#' and meta coordinates out side the circle
#' 
#' @param coord data.frame,spatial coordinates, 
#' should have column names "x" and "y" and "barcode".
#' @param km_coord_list The meta coordinates list. 
#' Generated by \code{\link{generate_kmeans_coord}}.
#' @param center the barcode of the center spot. 
#' @param radius_square radius of the circle. Default radius is 500 (micrometer).
#' The square is used to reduce steps and complexity of calculation. 
#' So the actual default value used in function is 500^2.
#'
#' @details This function will generate a new coordinates
#' which is a combination of meta and single spot coordinates.
#' In the circle of given center and radius, the coordinates will be keep 
#' as the orignal spot coordinates, and the outer will be replaced by the meta clustered coordinates.  
#' This step is to reduce the dimension of the data and fasten the calculation.
#' 
#' @return data.frame with coordinates and unique id.
#' @export 
generate_holed_coord <- function(
  coord,
  km_coord_list,
  center,
  radius_square=500^2
){
  meta_coord <- km_coord_list[[1]]
  km_clu_info  <- km_coord_list[[2]]
  coord_mat <- coord[,c("x","y")] %>% as.matrix()

  spatial_sub_mat <- 
    rep(1,nrow(coord_mat)) %x% coord_mat[center,,drop=FALSE] - coord_mat
  coord$dist_sq <- spatial_sub_mat^2 %>% rowSums()
  coord$position <- ifelse(coord$dist_sq <= radius_square, "inner", "outer")

  temp <- which(coord$position=="inner")
  holed_clu <- km_clu_info$cluster[temp] %>% unique()
  inner_clu_idx <- which(km_clu_info$cluster %in% holed_clu)

  inner_part <- coord[inner_clu_idx,]
  #inner_part$uni_id <- inner_part$barcode
  inner_part <- inner_part[,c("x","y","barcode")]

  outer_idx <- which(!(meta_coord$barcode %in% holed_clu))
  outer_part <- meta_coord[outer_idx,]
  #outer_part$uni_id <- outer_part$cluster
  outer_part <- outer_part[,c("x","y","barcode")]

  final_coord <- 
    rbind.data.frame(inner_part,outer_part)

  return(final_coord)
}

#' The implation of Generate Meta Expression
#'
#' The implation of generating a meta expression data.frame of clustered meta coordinates 
#' 
#' @param expr expression matrix or data.frame, commonly the result of generate_complex_data
#' @param clu_info cluster info, must be contain "barcode" and "cluster". 
#' The cluster is the result of Kmenas. 
#'
#' @details This function will sum-up the expresion of spot in a cluster. 
#' For receptor complex, we firstly calculate the complex expression of each spot by
#' generate_complex_data before this function.
#' @import dplyr
#' 
#' @return data.frame, sum-uped meta expression.
#' @export 
generate_meta_expr_impl <- function(
    expr,
    clu_info
){ 
  expr$barcode <- rownames(expr)
  expr %<>% left_join(clu_info[,c("barcode","cluster")],by="barcode")
  meta_expr <- expr %>%
    dplyr::group_by(cluster) %>%
    dplyr::select(!barcode) %>%
    dplyr::summarise(
      dplyr::across(everything(),sum)
    ) %>% 
    as.data.frame()
  rownames(meta_expr) <- meta_expr$cluster
  meta_expr$cluster <- NULL
  return(meta_expr)
}

#' Generate Holed Expression with coordinates
#'
#' Generate a holed expression with coordinates by concatenating 
#' holed coordinates and expression for given center and radius
#' 
#' @param all_expr all expression including each spot and meta expression.
#' @param radius radius of the circle. Default is 500 (micrometer).
#' In the circle of given center and radius, the coordinates will be kept 
#' as the orignal spot coordinates, and the outer will be replaced by the meta clustered coordinates.
#' @inheritParams generate_holed_coord
#' 
#' @seealso \code{\link{generate_holed_coord}}
#' 
#' @return data.frame, holed expression
#' @export 
generate_holed_coord_expr <- function(
  coord,
  km_coord_list,
  all_expr,
  center,
  radius=500
){
  holed_coord <- 
    generate_holed_coord(
      coord=coord,
      km_coord_list = km_coord_list,
      center=center,
      radius_square=radius^2
    )
  temp <- which(rownames(all_expr) %in% holed_coord$barcode)
  holed_expr <- cbind.data.frame(holed_coord,all_expr[temp,])
  return(holed_expr)
}

#'lapply warpper 
#' 
#' A warpper of lapply series functions, automatically select the proper functions.
#' 
#' @param seq_obj sequence object to be processed, such as list or vector
#' @param func function to be applied
#' @param verbose logical, whether to print progress, default is TRUE.
#' @param backend character, must be one of "auto","none","future", 
#' the backend of parallel computing, default is "auto".
#' @param ... arguments to be passed to lapply functions
#' 
#' @seealso \code{\link{lapply}} 
#' @seealso \code{\link[pbapply]{pblapply}}
#' @seealso  \code{\link[future.apply]{future_lapply}}
#' 
#' @import pbapply
#' @import future
#' @import future.apply
auto_select_lapply <- function(seq_obj,func,verbose=TRUE,backend=c("auto","none","future"),...){
  if(!is.vector(seq_obj)){stop("Must give a vector-like object to be processed")}
  para_set <- class(future::plan())
  backend <- match.arg(backend)

  if(backend=="auto"){
    if(!("sequential" %in% para_set) && ("multiprocess" %in% para_set)){
      backend <- "future"
    }else{
      backend <- NULL
    }
  }else if(backend=="none"){
      backend <- NULL
  }
  
  if(verbose){
    list_result <-pbapply::pblapply(X=seq_obj,FUN = func,cl=backend,...)
  }else if(is.null(backend)){
    list_result <-lapply(X=seq_obj,FUN = func,...)
  }else{
    list_result <-future.apply::future_lapply(X=seq_obj,FUN = func,...)
  }
  return(list_result)
}

# ' Simple Replicate one-col dgCMatrix
# ' 
# ' Simple Replicate one-col dgCMatrix to n-col
# ' 
# ' @param sparse_col one-col dgCMatrix
# ' @param n total number to be extend
# ' 
# ' 
# ' @importFrom Matrix sparseMatrix
# ' @importClassesFrom Matrix dgCMatrix CsparseMatrix
# ' @importFrom methods new
# replicate_sparse_column <- function(sparse_col, n) {
#   # Check if input is a single-column sparse matrix
#   if (!inherits(sparse_col, "dgCMatrix") || ncol(sparse_col) != 1) {
#     stop("Input must be a single-column dgCMatrix.")
#   }
  
#   # Create the replicated sparse matrix
#   # replicated <- Matrix::sparseMatrix(
#   #   i = rep(sparse_col@i, n),           # Replicate row indices
#   #   p = seq(0, length(sparse_col@x) * n, length.out = (n + 1)),      # Assign columns sequentially
#   #   x = rep(sparse_col@x, n),                   # Replicate the non-zero values
#   #   index1 = FALSE,  #fit the directly extraction form dgc
#   #   dims = c(nrow(sparse_col), n), # Set dimensions
#   #   repr = "C"  #dgC
#   # )
#   replicated <- new("dgCMatrix",
#                 i = rep(sparse_col@i, n) %>% as.integer(),               # Replicate row indices
#                 p = seq(0, length(sparse_col@x) * n, length.out = (n + 1)) %>% as.integer(), # Compute column pointers
#                 x = rep(sparse_col@x, n),               # Replicate non-zero values
#                 Dim = c(nrow(sparse_col), n) %>% as.integer()            # Dimensions
#   )
#   return(replicated)
# }

#' compute LR kept Matrix
#' 
#' Compute LR kept Matrix from one-col dgCMatrix and return dgCMatrix
#' 
#' @param q_net_sign one-col dgCMatrix
#' 
#' @importFrom Matrix sparseMatrix
#' @importClassesFrom Matrix dgCMatrix CsparseMatrix
#' 
#' @return a dgCMatrix store the LR kept matrix
compute_LR_kept_dgC <- function(q_net_sign) {
  # Check if input is a single-column dgCMatrix
  if (!inherits(q_net_sign, "dgCMatrix") || ncol(q_net_sign) != 1) {
    stop("Input must be a single-column dgCMatrix.")
  }
  
  # Extract non-zero elements and their row indices
  x <- q_net_sign@x
  i <- q_net_sign@i
  n <- nrow(q_net_sign)
  
  # Generate all pairwise differences as a triplet representation
  row_idx <- rep(i, each = length(x))    # Row indices
  col_idx <- rep(i, times = length(x))  # Column indices
  values <- rep(x, each = length(x)) - rep(x, times = length(x)) # Pairwise differences
  
  values[values<2] <- 0
  values[values==2] <- 1

  # Filter out zeros
  nonzero <- values != 0
  row_idx <- row_idx[nonzero]
  col_idx <- col_idx[nonzero]
  values <- values[nonzero]
  
  # Construct the result as a dgCMatrix
  result <- Matrix::sparseMatrix(,
    i = row_idx %>% as.integer(),
    j = col_idx %>% as.integer(),
    x = values,
    dims = c(n, n),
    index1 = FALSE,
    giveCsparse = TRUE                            # Matrix dimensions
  ) 
  return(result)
}
