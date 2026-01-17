#**********************************

#**************************
#** Tensor Decomposition **
#**************************


#' prepare tensor data for CP decomposition
#'
#' This function pre-processes the LR pair field data from a spWAVE object
#'
#' @param spwave spWAVE object containing LR pair field data
#' @param low_quantile quantile threshold to filter low expression values. Default is 0.1 (10%).
#'
#' @return log-normalized filtered tensor array for CP decomposition
#' @export
prep_TD_data <- function(spwave,low_quantile=0.1){
  # create tensor data
  tensor_data <- spwave@LR_pair_field

  tensor_data_log <- array(data=NA,
    dim=c(dim(tensor_data)[1],3,dim(tensor_data)[3]),
    dimnames = list(dimnames(tensor_data)[[1]],c("logNorm","logEx","logEy"),dimnames(tensor_data)[[3]])
  )

  tensor_data_log[,"logNorm",] <-
    log(tensor_data[,"Ex",]^2+tensor_data[,"Ey",]^2)/2 #log 歪除2等于取了开方，节省计算开销
  tensor_data_log[,"logEx",] <-
    abs(tensor_data_log[,"logNorm",])*cos(atan2(tensor_data[,"Ey",],tensor_data[,"Ex",]))
  tensor_data_log[,"logEy",] <-
    abs(tensor_data_log[,"logNorm",])*sin(atan2(tensor_data[,"Ey",],tensor_data[,"Ex",]))

  # filter low expression
  tensor_norm <-
    apply(tensor_data,3,function(x){
      sqrt(x[,"Ex"]^2+x[,"Ey"]^2)
    })

  low_cutoff <- tensor_norm %>% as.vector %>% quantile(probs=c(low_quantile))

  cond <- tensor_data_log[, "logNorm", ] < log(low_cutoff)
  cond_expanded <- array(cond, dim = c(dim(cond), ncol(tensor_data_log)))
  cond_expanded <- aperm(cond_expanded, perm = c(1, 3, 2))
  tensor_data_log[cond_expanded] <- 0
  return(tensor_data_log[,2:3,])
}


#' Repeated Tensor CP Decomposition
#'
#' This function performs CP decomposition on a tensor multiple times for a range of ranks
#'
#' @param tensor Input tensor array for CP decomposition
#' @param num_rep Number of repetitions for each rank. Default is 10.
#' @param rank_range Range of ranks to test for CP decomposition. Default is 10:15.
#' @param random_seed Seed for random number generation to ensure reproducibility. Default is 42.
#' @param const character vextor, Constraints for each mode during CP decomposition. Default is c("uncons","uncons","orthog").
#' see multiway::parafac and CMLS::cmls for more details and available options.
#' @param corcondia_cutoff Minimum corcondia value to accept a decomposition. Default is 10.
#' @param attempts Maximum number of attempts to find a valid decomposition. 
#' Default is 7 and must not be lower than 3.
#' @param parallel Logical indicating whether to use parallel processing. Default is TRUE.
#' @param cl Cluster object for parallel processing. Default is NULL.
#'
#'
#' @seealso \code{\link[multiway]{parafac}}
#' @seealso \code{\link[CMLS]{cmls}}
#' @seealso \code{\link[multiway]{corcondia}}
#'
#' @return nested list of CP decomposition results for each rank
#' @export
repeat_tensor_CPD <- function(
  tensor,
  num_rep=10,
  rank_range=10:15,
  random_seed=42,
  const=c("uncons","uncons","orthog"),
  corcondia_cutoff=10,
  attempts=7,
  parallel=TRUE,cl=NULL){
  set.seed(random_seed)

  if (!requireNamespace("multiway", quietly = TRUE)) {
    stop(
      "Package \"multiway\" must be installed to use this function.",
      call. = FALSE
    )
  }


  if(parallel){
    if (!requireNamespace("parallel", quietly = TRUE)) {
    stop(
      "Package \"parallel\" must be installed to use this function.",
      call. = FALSE
    )
  }
    if(!is.null(cl)){
      parallel::clusterEvalQ(cl,library(multiway))
    }else{
      stop("Please provide a valid cluster object 'cl' for parallel processing.")
    }
  }

  CP_decomp <- list()
  loop_idx <- seq_along(rank_range)
  lowering_label <- FALSE
  for(i in loop_idx){
    R=rank_range[i]


    if(i!=1 && sum(sapply(CP_decomp[[i-1]], is.null)) >3) {
      cat("Stop decompostion after rank", R,  "due to invaild result.\n")
      break
    }

    if(lowering_label){
      rank_cutoff <- current_cutoff + 1
    }else{
      rank_cutoff <- corcondia_cutoff
    }
    cat("Processing rank", R, ":\n")
    current_cutoff <- rank_cutoff
    rep_results <- list()

    for(j in 1:num_rep) {
      cat("Repeat", j, "...\n")

      reject_counter <- 0
      #for(atmpts in 1:attempts){
      while(reject_counter < attempts){
        temp_result <- multiway::parafac(tensor,nfac=R,
            const=const,
            parallel=parallel,cl=cl)
        corcondia_temp <- multiway::corcondia(tensor, temp_result)
        print(corcondia_temp)

        if(corcondia_temp > current_cutoff){
          rep_results[[j]] <- temp_result
          reject_counter <- 0
          break
        }else{
          reject_counter <- reject_counter + 1
          cat("corcondia below", current_cutoff, "try again...\n")
          
          if(reject_counter >= 3){
            if(current_cutoff > 1){
              if(corcondia_temp > 2){
                current_cutoff <- max((current_cutoff + corcondia_temp)/2,1)
              }else{
                current_cutoff <- max(current_cutoff - 1,1)
              }
              
              reject_counter <- 0
              lowering_label <- TRUE
              cat("Lowering corcondia threshold to", current_cutoff, "try again...\n")
            }else{
              rep_results[[j]] <- NULL
              cat("No valid result, return NULL\n")
              break
            }

          }
        }
      }
    }
    CP_decomp[[i]] <- rep_results
  }
  names(CP_decomp) <- paste0("Rank_", rank_range[1:length(CP_decomp)])
  return(CP_decomp)
}



#' Construct Correlation Matrix from Tensor Decomposition Results
#'
#' This function constructs a correlation matrix from tensor decomposition results
#' by calculating pairwise correlations between components and applying significance
#' and effect size thresholds.
#'
#' @param TD_result Tensor decomposition results from CP decomposition. Can be either
#'   a nested list structure or a direct list of decomposition results.
#' @param R_cutoff Numeric value specifying the minimum absolute correlation coefficient
#'   to retain. Correlations with absolute value below this threshold will be set to 0.
#'   Default is 0.5.
#' @param p_cutoff Numeric value specifying the maximum p-value for significance.
#'   Correlations with p-values above this threshold will be set to 0 after FDR correction.
#'   Default is 0.05.
#' @param link_cutoff Numeric value specifying the minimum number of effective correlations with other components. 
#'   the components with less effective correlations lower than this value will be dropped.
#'   Default is 10.
#' @param return_raw Logical indicating whether to return additional information beyond
#'   the filtered correlation matrix. If TRUE, returns a list containing the filtered
#'   correlation matrix, raw correlation matrix, and adjusted p-values. 
#'   If FALSE, returns only the filtered correlation matrix. Default is TRUE.
#'
#' @return If \code{return_raw = FALSE}, returns a numeric matrix of filtered
#'   correlations. If \code{return_raw = TRUE}, returns a list containing:
#'   \item{cor_ft}{Filtered correlation matrix with non-significant correlations set to 0}
#'   \item{cor_raw}{Raw correlation matrix before filtering}
#'   \item{p_adj}{Matrix of FDR-adjusted p-values}
#'
#' @examples
#' \dontrun{
#' # Example with tensor decomposition results
#' # Assuming TD_result contains CP decomposition results
#' cor_matrix <- construct_cor_mat(TD_result, R_cutoff = 0.6, p_cutoff = 0.01)
#'
#' # Return full results including raw correlations and p-values
#' full_results <- construct_cor_mat(TD_result, return_raw = TRUE)
#' str(full_results)
#' }
#'
#' @export
construct_cor_mat <- function(
  TD_result,
  R_cutoff=0.5,
  p_cutoff=0.05,
  link_cutoff=10,
  return_raw=TRUE){

  if (!requireNamespace("Hmisc", quietly = TRUE)) {
    stop(
      "Package \"Hmisc\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(is.null(TD_result[[1]]$C)){
    CP_decomp <- unlist(TD_result,recursive = FALSE)
  }else{
    CP_decomp <- TD_result
  }

  for(r in seq_along(CP_decomp)){
    H <- CP_decomp[[r]]$C
    colnames(H) <- paste0("R",ncol(H),"rep",r,"_",1:ncol(H))
    CP_decomp[[r]]$MP_temp <- as.data.frame(H)
  }

  MP_temp_list <- lapply(CP_decomp,function(x) x$MP_temp)
  MP_cbind_temp <- do.call(cbind,MP_temp_list)

  rc <- Hmisc::rcorr(as.matrix(MP_cbind_temp), type = "pearson")
  p_adj  <- matrix(p.adjust(rc$P, "fdr"), nrow = nrow(rc$P))
  sig <- which(!(p_adj < p_cutoff & abs(rc$r) > R_cutoff), arr.ind = TRUE)
  r_ft <- rc$r
  r_ft[sig] <- 0

  r_ft_bin <- r_ft
  r_ft_bin[r_ft_bin != 0] <- 1
  temp <- which(rowSums(r_ft_bin) <= link_cutoff)
  r_ft <- r_ft[-temp,-temp]

  if(return_raw){
    return(list(cor_ft=r_ft,cor_raw=rc$r,p_adj=p_adj))
  }else(
    return(r_ft)
  )
}





#' Dynamic Module Assignment from Hierarchical Clustering
#'
#' Dynamic module assignment based on hierarchical clustering. This function assigns
#' modules to leaves based on their distance from the root of the dendrogram and
#' ensures large clusters (≥MIN_THRESHOLD) are not merged into larger clusters.
#'
#' @param hc hclust object representing the hierarchical clustering
#' @param offset distance offset from the root to control cutting height. Default is NULL and will be set to max height.
#' @param MIN_THRESHOLD the minimum size of clusters to prevent merging, 
#'        smaller threshold will create smaller size of cluster, more unassign leaves, and more clusters. 
#'        Default is 10.
#' @param Z_THRESHOLD distance threshold, controls the strictness of cutting, 
#'        smaller threshold will create larger size of cluster, less unassign leaves, and more clusters. 
#'        Default is the median of hc$height.
#'
#' @return named integer vector mapping leaf node labels to cluster IDs (-1 indicates unassigned)
#'
#' @examples
#' \dontrun{
#' 
#' hc <- hclust(dist_mat, method = "average")
#' clusters <- adynamic_module_cutree(hc)
#'
#' clusters <- adynamic_module_cutree(hc,MIN_THRESHOLD = 20, Z_THRESHOLD = 5)
#' }
#' @export
dynamic_module_cutree <- function(hc,
  offset = NULL,
  MIN_THRESHOLD = 10,
  Z_THRESHOLD = median(hc$height)) {
  if (!inherits(hc, "hclust")) {
    stop("hc must be an hclust object")
  }

  n_leaves <- length(hc$order)
  n_nodes <- nrow(hc$merge)

  if (is.null(offset)) {
    offset <- max(hc$height)
  }

  # Get each node's member count
  node_memberships <- compute_node_memberships(hc$merge)

  labels <- rep(-1, n_nodes)
  clust_counter <- 0  

  for (i in seq_len(n_nodes)) {
    left_child <- hc$merge[i, 1]
    right_child <- hc$merge[i, 2]

    # Get left subtree member count and current label
    if (left_child < 0) {
      left_size <- 1
      left_label <- -1
    } else {
      left_size <- node_memberships[abs(left_child)]
      left_label <- labels[abs(left_child)]
    }

    # Get right subtree member count and current label
    if (right_child < 0) {
      right_size <- 1
      right_label <- -1
    } else {
      right_size <- node_memberships[abs(right_child)]
      right_label <- labels[abs(right_child)]
    }


    if (left_size >= MIN_THRESHOLD && right_size >= MIN_THRESHOLD) {
      # subtrees too large，do not merge
      new_label <- -1
    } else if (hc$height[i] > (offset - Z_THRESHOLD)) {
      # distance too large, do not merge
      new_label <- -1
    } else if (left_size >= MIN_THRESHOLD) {
      # merge into left subtree
      new_label <- left_label
    } else if (right_size >= MIN_THRESHOLD) {
      # merge into right subtree
      new_label <- right_label
    } else if ((left_size + right_size) >= MIN_THRESHOLD) {
      # merge and form new cluster
      new_label <- clust_counter
      clust_counter <- clust_counter + 1
    } else {
      new_label <- -1
    }

    labels[i] <- new_label
  }

  out_clusters <- rep(-2, n_leaves)

  # spread labels from root to leaves
  if (n_nodes > 0) {
    root_label <- labels[n_nodes]
    out_clusters <- prop_label(n_nodes, root_label, labels, out_clusters, hc$merge)
  } else {
    out_clusters[1] <- 0
  }

  # safety check
  unassigned <- out_clusters == -2
  if (any(unassigned)) {
    warning(sprintf("%d leaves were not assigned to any cluster", sum(unassigned)))
    out_clusters[unassigned] <- -1
  }

  unique_clusters <- sort(unique(out_clusters))
  unique_clusters <- unique_clusters[unique_clusters != -1]

  if (length(unique_clusters) > 0) {
    clust_map <- stats::setNames(seq_len(length(unique_clusters)), unique_clusters)
    clust_map <- c(clust_map, "-1" = -1)
    out_clusters <- clust_map[as.character(out_clusters)]
  }

  names(out_clusters) <- hc$labels

  return(out_clusters)
}

#' Compute Node Memberships Number in Hierarchical Clustering
#'
#' @param merge_mat hclust object's merge matrix, e.g. hc$merge
#' @return integer vector of node memberships count
#'
compute_node_memberships <- function(merge_mat) {
  n_nodes <- nrow(merge_mat)
  memberships <- integer(n_nodes)

  for (i in seq_len(n_nodes)) {
    left_child <- merge_mat[i, 1]
    right_child <- merge_mat[i, 2]

    left_size <- if (left_child < 0) {
      1
    } else {
      memberships[abs(left_child)]
    }

    right_size <- if (right_child < 0) {
      1
    } else {
      memberships[abs(right_child)]
    }

    memberships[i] <- left_size + right_size
  }

  return(memberships)
}

#' Spread Cluster Labels from Internal Nodes to Leaf Nodes
#'
#' @param node_id the node ID to start propagation (1-based index for internal nodes)
#' @param label the label to inherit (-1 means use the node's own label)
#' @param labels internal node labels vector
#' @param out_clusters leaf node output vector (indexed by leaf node IDs)
#' @param merge_mat hclust object's merge matrix
#'
#' @return updated out_clusters with new labels
#'
prop_label <- function(node_id, label, labels, out_clusters, merge_mat) {
  current_label <- if (label == -1) labels[node_id] else label

  left_child <- merge_mat[node_id, 1]
  right_child <- merge_mat[node_id, 2]

  if (left_child < 0) {
    leaf_idx <- abs(left_child)
    out_clusters[leaf_idx] <- current_label
  } else {
    out_clusters <- prop_label(abs(left_child), current_label,
                                labels, out_clusters, merge_mat)
  }

  if (right_child < 0) {
    leaf_idx <- abs(right_child)
    out_clusters[leaf_idx] <- current_label
  } else {
    out_clusters <- prop_label(abs(right_child), current_label,
                                labels, out_clusters, merge_mat)
  }

  return(out_clusters)
}


#' Cluster Modules from Correlation Matrix
#' 
#' This function performs hierarchical clustering on a correlation matrix and assigns modules
#' based on specified thresholds.
#' 
#' @param cor_mat Correlation matrix to cluster
#' @param linkage_method Linkage method for hierarchical clustering. Default is "average". 
#'        Check stats::hclust for more options.
#' @param Z_THRESHOLD Distance threshold to control the strictness of cutting. 
#'        Default is NULL and will be set to median of hc$height.
#' @inheritParams dynamic_module_cutree
#' 
#' @seealso \code{\link[stats]{hclust}}
#' @seealso \code{\link{dynamic_module_cutree}}
#' 
#' @return A list containing:
#' \item{hc}{hclust object representing the hierarchical clustering}
#' \item{module}{characters vector mapping leaf node labels to cluster IDs}
#' \item{module_df}{Two cols data.frame with module assignments for each leaf node}
#' @export
#' 
perform_dynamic_module_clustering <- function(
  cor_mat,
  linkage_method="average",
  offset=NULL,
  MIN_THRESHOLD=10,
  Z_THRESHOLD=NULL

){
  dist_matrix <- stats::as.dist(1 - cor_mat)
  hc <- stats::hclust(dist_matrix, method = linkage_method)

  if(is.null(Z_THRESHOLD)){
    Z_THRESHOLD <- median(hc$height)
  }

  module_assign <- dynamic_module_cutree(
    hc,
    offset=offset,
    MIN_THRESHOLD=MIN_THRESHOLD,
    Z_THRESHOLD=Z_THRESHOLD
  )

  module_df <- as.data.frame(module_assign)
  colnames(module_df)[1] <- "module"
  module_df$module <- paste0("module_",module_df$module)
  module_df$member <- rownames(module_df)

  module <- split(module_df[,"member"],module_df$module)
  names(module)[1] <- "other"
  module_df[module[[1]],"module"] <- "other"

  return(list(
    hc=hc,
    module=module,
    module_df=module_df
  ))
}

#' Plot Correlation Matrix with Module Annotations
#' 
#' @param cor_mat Correlation matrix to plot  
#' @param module_list Output of clustering_modules function
#' @param ht_col Color palette for heatmap. Default is NULL 
#'        and will be set to colorRamp2(c(-1, 0, 1), c("navy", "white", "firebrick3")).
#' @param module_col Color palette for module annotations. Default is NULL 
#'        and will be set to MD2_color_picker(length(module_list$module)-1).
#' @param ... Additional arguments to pass to ComplexHeatmap::Heatmap
#' 
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' 
#' @export
plot_cor_module_heatmap <- function(
  cor_mat,
  module_list,
  ht_col=NULL,
  module_col=NULL,
  ...
){
  #prepare ht data
  hc_order <- module_list$hc$labels
  module_annotation <- module_list$module_df

  draw_ht_data <- cor_mat[hc_order, hc_order]

  #prepare colors
  if(is.null(ht_col)){
    ht_col <- circlize::colorRamp2(
      breaks = c(-1, 0, 1),
      colors = c("navy", "white", "firebrick3")
    )
  }#else #*TODO 后续需要添加 col_fun checker


  if(is.null(module_col)){
    module_col <- MD2_color_picker(length(module_list$module)-1)
    names(module_col) <- names(module_list$module)[2:length(module_list$module)]
    module_col["other"] <- "#D9D9D9"
  }#else #*TODO 后续需要添加 module_col checker， 检查长度和names

  #prepare module annotation
  module_ha <- ComplexHeatmap::HeatmapAnnotation(
    Module = factor(module_annotation$module,levels=names(module_col)),
    col = list(Module=module_col),
    annotation_name_side = "left",
    show_legend = TRUE,
    annotation_legend_param = list(
      Module = list(
        title = "Modules",
        ncol = 1,
        labels = names(module_col)
      )
    )
  )

  ComplexHeatmap::Heatmap(
    matrix = draw_ht_data,
    col = ht_col,
    cluster_rows = module_list$hc,
    cluster_columns = module_list$hc,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    top_annotation = module_ha,
    name = "Pearson\nCorrelation",
    use_raster=TRUE,
    column_title = "Hierarchical Clustering with Module Annotation",
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1)  # 自定义图例刻度
    ),
    ...
  )

}


#' Extract Top Ranked Members from Modules
#' 
#' This function extracts the top ranked members from each module based on their
#' contributions in the tensor decomposition results.
#' 
#' @param TD_result Tensor decomposition results from CP decomposition. Can be either
#'   a nested list structure or a direct list of decomposition results.
#' @param tensor filtered tensor array used for CP decomposition, generated by prep_TD_data()
#' @param module_list Output of clustering_modules function
#' @param top_n Number of top-ranked members to extract from each module. Default is 50.
#' 
#' @return A list of top-ranked members for each module
#' @export
extract_module_topLR <- function(
  TD_result,
  tensor,
  module_list,
  top_n=50
){
  if(is.null(TD_result[[1]]$C)){
    CP_decomp <- unlist(TD_result,recursive = FALSE)
  }else{
    CP_decomp <- TD_result
  }

  for(r in seq_along(CP_decomp)){
    H <- CP_decomp[[r]]$C
    rownames(H) <- dimnames(tensor)[[3]]
    colnames(H) <- paste0("R",ncol(H),"rep",r,"_",1:ncol(H))
    CP_decomp[[r]]$MP_temp <- as.data.frame(H)
  }

  #Normalize and ignore the negative values
  CP_decomp %<>% lapply(function(x){
    x$MP_norm <- apply(x$MP_temp,2,function(y){
      value=y/sqrt(sum(y^2))
      value[which(value<=0)] <- NA
      return(value)
    })
    return(x)
  })

  top_list <- lapply(CP_decomp,function(x){
    apply(x$MP_norm,2,function(y){
      names(sort(y,decreasing = T)[1:top_n])
    }) %>% as.data.frame()
  })

  top_flat_df <- do.call(cbind,top_list)

  module_member <- list()
  for(i in 2:length(module_list$module)){ #* the first is "other"
    module_name <- names(module_list$module)[i]
    member_vec <- module_list$module[[i]]
    member_indices <- match(member_vec,colnames(top_flat_df))
    sub_top_df <- top_flat_df[,member_indices]
    sub_top_vec <- as.vector(as.matrix(sub_top_df))

    # calc rank scores for genes in each component, then sum-up score
    rank_scores <- list()
    for(col in 1:ncol(sub_top_df)){
      for(row in 1:nrow(sub_top_df)){
        gene <- sub_top_df[row, col]
        if(!is.na(gene)){
          rank_scores[[gene]] <- rank_scores[[gene]] %||% 0
          rank_scores[[gene]] <- rank_scores[[gene]] + (top_n+1 - row)  # higher rank gets more points
        }
      }
    }

    # Rank sum-up scores and get top n
    rank_scores_vec <- unlist(rank_scores)
    sorted_genes <- sort(rank_scores_vec, decreasing = TRUE)

    final_members <- names(sorted_genes)[1:min(top_n, length(sorted_genes))]

    module_member[[module_name]] <- final_members
  }

  return(module_member)
}


#' Scoring Cell Modules using ssGSEA
#' 
#' This function scores cell modules using single-sample Gene Set Enrichment Analysis (ssGSEA)
#' 
#' @param tensor filtered tensor array used for CP decomposition, generated by prep_TD_data()
#' @param module_LR_list A list of LR pairs for each module. Generated by extract_module_topLR()
#' @param ... Additional arguments to pass to GSVA::gsva, e.g. BPPARAM, verbose, etc.
#' 
#' @details ssgesa scoring is performed using the GSVA package. The kcdf is set to "Gaussian", 
#'          since the expression matrix is log-normalized. Other parameters can be adjusted.
#' 
#' @return A matrix of ssGSEA scores for each module across samples
#' @examples 
#' \dontrun{
#' 
#' # Assuming tensor and module_LR_list are prepared
#' ssGSEA_scores <- scoring_cell_module(tensor, module_LR_list)
#' 
#' # Using additional parameters such as parallel.sz, BPPARAM, and verbose for GSVA 
#' library(BiocParallel)
#' ssGSEA_scores <- scoring_cell_module(
#'   tensor,
#'   module_LR_list,
#'   parallel.sz = 16,
#'   BPPARAM = MulticoreParam(workers = 16,tasks=16, progressbar = T), 
#'   # use SnowParam for Windows system.
#'   verbose = TRUE
#' )
#' 
#' }
#' @export

scoring_cell_module <- function(
  tensor,
  module_LR_list,
  ...
){
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop(
      "Package \"GSVA\" must be installed to use this function.",
      call. = FALSE
    )
  }
  # prepare expression matrix
  expr_mat <- apply(tensor,3,function(x){
    sqrt(x[,"Ex"]^2+x[,"Ey"]^2)
  })
  rownames(expr_mat) <- dimnames(tensor)[[1]]
  colnames(expr_mat) <- dimnames(tensor)[[3]]

  # run ssGSEA
  ssGSEA_result <- GSVA::gsva(
    expr = expr_mat,
    gset.idx.list = module_LR_list,
    method = "ssgsea",
    kcdf = "Gaussian",
    ...
  )

  return(ssGSEA_result)
}


#' Assign Cell Modules based on ssGSEA Scores
#' 
#' This function assigns cell modules based on ssGSEA scores with specified cutoffs and
#' switch the assignments for cells in small clusters or with low score differences.
#' 
#' @param ssGSEA_result A matrix of ssGSEA scores for each module across samples
#' @param score_cutoff Minimum ssGSEA score to assign a module. Default is 0.1.
#' @param score_diff_cutoff Minimum score difference between top two modules to keep assignment. Default is 0.03.
#'        Set 0 to keep all assignments regardless of minimum cell numbers of clusters.
#' @param min_cells Minimum number of cells required for a module to be retained. Default is 0.5% of total cells.
#' 
#' @return A data frame with cell barcodes and their assigned modules
#' @export
assign_cell_module <- function(
  ssGSEA_result,
  score_cutoff=0.1,
  score_diff_cutoff=0.03,
  min_cells=0.005*ncol(ssGSEA_result)
){
  assign_modules <- apply(ssGSEA_result,2,function(x){
    if(max(x)>score_cutoff){
      assign <- names(which(x == max(x)))
    }else{
      assign <- "other"
    }
    return(assign)
  })

  # If score_diff_cutoff is 0, return assignments directly
  if(score_diff_cutoff == 0){
    final_assign <- data.frame(
      barcode=names(assign_modules),
      assign=assign_modules,
      row.names = names(assign_modules)
    )
    return(final_assign)
  }

  #If score_diff_cutoff > 0, switch assignments for small clusters with low score diff
  assign_modules_2nd <- apply(ssGSEA_result,2,function(x){
    if(max(x)>0.1){
      assign <- names(x)[order(x,decreasing = T)[2]]
    }else{
      assign <- "other"
    }
    return(assign)
  })

  assign_diff <- apply(ssGSEA_result,2,function(x){
      sorted_x <- x[order(x,decreasing = T)]
      diff <- sorted_x[1] - sorted_x[2]
    return(diff)
  })

  assign_df <-
    cbind.data.frame(assign_modules,assign_modules_2nd,assign_diff)

  temp <- table(assign_df[,1])
  dropped_cand <- names(temp)[temp<min_cells]
  dropped_cell <-
    rownames(assign_df)[which(assign_df[,1] %in% dropped_cand)]

  switched <- sapply(dropped_cell,function(id){
    if(as.numeric(assign_df[id,3]) < score_diff_cutoff &&
      !(assign_df[id,2] %in% dropped_cand) ){
        return(assign_df[id,2])
    }else{
      return("other")
    }
  },simplify = T)

  assign_df$assign <- assign_df$assign_modules
  assign_df[dropped_cell,"assign"] <- switched

  final_assign <- data.frame(
    barcode=rownames(assign_df),
    assign=assign_df$assign,
    row.names = rownames(assign_df)
  )

  return(final_assign)
}
