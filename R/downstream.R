#*********************
#** tradeSeq Module **
#*********************

#' Prepare data for tradeSeq analysis
#'
#' Extract and filter expression and field result for downstream tradeSeq analysis
#'
#' @param seurat_obj Seurat Object.
#' @param field_df a data.frame, the field result of interest.
#' @param traj_arg character, must be a colname in field result data.frame.
#' @param min_expr minimum expression in counts. default is 1
#' @param min_n_cell minimum number of cell expressing a gene. default NULL. This arg will mask 'min_pct_cell'.
#' @param min_pct_cell minimum percentage of cell expressing a gene. Default is 0.01
#'
#' @return return a list including filtered expression matrix in counts, 
#' a data.frame contain barcode and weight which used in tradeSeq.
#' @import Seurat
#' @export
prep_tradeSeq_data <- function(seurat_obj,field_df,traj_arg,min_expr=1,
  min_n_cell=NULL,
  min_pct_cell=0.01){
  cell_PT_df <- field_df[,traj_arg,F]
  cell_PT_df$cellWeight <- 1

  expr_mat <- GetAssayData(seurat_obj,layer="counts")
  expr_mat_ft <- filter_expr_by_cutoff(expr_mat,
    min_expr=min_expr,
    min_n_cell=min_n_cell,
    min_pct_cell=min_pct_cell)
  
  return(list(expr_mat=expr_mat_ft,cell_PT=cell_PT_df))
}



#' Run TradeSeq for field trajectory DE analysis
#'
#' Run tradeSeq based on given field result, and calc the field trajectory relative DE.
#'
#' @param seurat_obj Seurat Object.
#' @param database_result field result generated from upstream.
#' @param LR LR pair or family you interest.
#' @param kept_db kept_db matching the database result.
#' @param field_df a data.frame, the field result of interest. Default is NULL.
#' This args will mask database_result, LR, and kept_db which generate field_df. See detail.
#' @param traj_arg character, must be a colname in field result data.frame.
#' @param min_expr minimum expression in counts. default is 1
#' @param min_n_cell minimum number of cell expressing a gene. default NULL. This arg will mask 'min_pct_cell'.
#' @param min_pct_cell minimum percentage of cell expressing a gene. Default is 0.01
#' @param nknots number of knots, passing to tradeSeq, Default is 6.
#' @param verbose logical, whether to show verbose.
#' @param parallel logical, wether to use parallel.
#' @param BPPARAM the BiocParallel program, Default is NULL, see detail.
#' @param paral_workers number of parallel workers. Default is 8. It depends on BiocParallel. 
#' You can set certain parallel program before this function. Please set same workers.
#'
#' @details This function will run tradeSeq from given Seurat object and field info to 
#' calc the field related different gene expression.  
#' The field info can generate directly from upstream result by the args \code{database_result},\code{LR}, and \code{kept_db}.
#' Or manually give a data.frame which you can add field related quantity and give the colnames of the quantity.
#' The data.frame must contain barcode, coordinates. \cr
#' The parallel program depends on the BiocParallel. 
#' The default program is different from platform, please check document of BiocParallel.
#' Or you can manually set the program before the function runs.
#' 
#' @return return tradeSeq result
#' @export
run_LR_tradeSeq <- function(
  seurat_obj,
  database_result=NULL,LR=NULL,kept_db=NULL,
  field_df=NULL,
  traj_arg="U",
  min_expr=1,min_n_cell=NULL,min_pct_cell=0.01,
  nknots=6,verbose=TRUE,parallel=TRUE,BPPARAM=NULL,paral_workers=8
){
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop(
      "Package \"tradeSeq\" must be installed to use this function.",
      call. = FALSE
    )
  }
  data_check <- 
    any(
      all(!is.null(database_result),!is.null(LR)),
      !is.null(field_df)
  )
  if(!data_check){
    stop("Must give database_result, LR, and kept_db, or manually give field_df.")
  }
  message("prep data")
  #** prep data
  if(is.null(field_df)){
    field_df <- 
      perform_field_extract(database_result=database_result,LR=LR,kept_db=kept_db)
    field_df %<>% calc_field_strength()
  }
  if(!(traj_arg %in% colnames(field_df))){
    stop(paste0(traj_arg," not found in field_df"))
  }
  trade_data <- 
    prep_tradeSeq_data(
      seurat_obj=seurat_obj,field_df=field_df,traj_arg = traj_arg,
      min_expr=min_expr,
      min_n_cell=min_n_cell,
      min_pct_cell=min_pct_cell)

  if(parallel && requireNamespace("BiocParallel", quietly = TRUE)){
    if(is.null(BPPARAM)){
      BPPARAM <- BiocParallel::bpparam()
    }else(BPPARAM = BPPARAM)
    BPPARAM$workers <- paral_workers
    print(BPPARAM)
  }else{
    BPPARAM=NULL
  }
  message("running tradeSeq")
  sce <- 
    tradeSeq::fitGAM(counts = trade_data$expr_mat, 
      pseudotime=trade_data$cell_PT[,traj_arg,F], 
      cellWeights=trade_data$cell_PT[,"cellWeight",F],
      nknots=nknots, 
      verbose=verbose,
      parallel=parallel, 
      BPPARAM = BPPARAM
      )

  return(list(trade_res=sce,field_df=field_df))
}


#' Run tradeSeq gene test
#'
#' Run tradeSeq gene test by given tradeSeq result and selected method
#'
#' @param trade_res tradeSeq result from \code{run_LR_tradeSeq}.
#' @param method_use tradeTest args，"association" or "startvsend", default is "association".
#'
#' @return return data.frame contain gene pvalue, p.adj, waldStat, without test
#' @details Please see tradeSeq documentation of \code{associationTest} and \code{startvsend} for more details.
#' @importFrom stats p.adjust
#' @import dplyr 
#' @export
run_tradeSeq_gene_test <- function(
    trade_res,
    method_use=c("association","startvsend")
){
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop(
      "Package \"tradeSeq\" must be installed to use this function.",
      call. = FALSE
    )
  }
  method_use=match.arg(method_use)
  tradeTest <- switch(method_use,
    association = tradeSeq::associationTest,
    startvsend = tradeSeq::startVsEndTest
  )
  testRes <- tradeTest(trade_res)
  testRes[is.na(testRes[,"waldStat"]),"waldStat"] <- 0.0
  testRes[is.na(testRes[,"df"]),"df"] <- 0.0
  testRes[is.na(testRes[,"pvalue"]),"pvalue"] <- 1.0
  testRes$p.adj <-stats::p.adjust(testRes$pvalue,method = 'BH')
  testRes %<>% dplyr::arrange(dplyr::desc(waldStat))
  return(testRes)
}



#' Plot gene expression in trajectory
#'
#' Plot heatmap of gene expression in trajectory, from tradeSeq result 
#' 
#' @param trade_res sce, result from \code{run_LR_tradeSeq}.
#' @param tradeTest_res data.frame,  result from \code{run_tradeSeq_gene_test}.
#' @param field_df field info given in the \code{run_LR_tradeSeq}.
#' @param p_val p value type to filter, "p.adj" or "pvalue", default is "p.adj".
#' @param p_cutoff p value cutoff, default is 0.05.
#' @param top_n_gene top n genes to show. Default is Inf.
#' @param order_by keyword to rank genes, "wald" or "logFC", default is "wald".
#' @param select_genes gene list to show, will mask top_n_gene and order_by.
#' @param title character, plot title.
#' @param ht_col heatmap color.
#' @param pot_col potential color, should be 100 breaks.
#' @param show_rownames logical, whether to show rownames. Default is FALSE.
#' @param n_split_row cut row into n_split_row tree. Default is 4.
#'
#' @return return plot
#' @import viridis
#' @import ComplexHeatmap
#' @import stringr
#' @import dplyr
#' @import grid
#' @importFrom circlize colorRamp2
#' @importFrom utils head
#' @export
plot_tradeTest_heatmap <- function(
  trade_res,
  tradeTest_res,
  field_df,
  p_val=c("p.adj","pvalue"),
  p_cutoff=0.05,
  top_n_gene=Inf,
  order_by=c("wald","logFC"),
  select_genes=NULL,
  title=NULL,
  ht_col=viridis::plasma(100),
  pot_col=NULL,
  show_rownames=FALSE,
  n_split_row=4
){
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop(
      "Package \"tradeSeq\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if(is.null(select_genes)){
    order_by <- match.arg(order_by)
    od_col <- ifelse(order_by=="wald",1,4)
    top_n_gene <- min(top_n_gene,nrow(tradeTest_res))
    p_val <- match.arg(p_val)
    select_genes <- 
      tradeTest_res %>% 
        dplyr::filter(.data[[p_val]] < p_cutoff) %>%
        dplyr::arrange(dplyr::desc(colnames(.data)[od_col])) %>% 
        head(top_n_gene) %>% 
        rownames()
  }else{
    if(!all(select_genes %in% rownames(tradeTest_res))){
      warning("Not all select_genes are in the result, using intersected instead.")
    }
    select_genes <- intersect(select_genes,rownames(tradeTest_res))
  }
  print(select_genes)
  yhatSmooth <- 
    tradeSeq::predictSmooth(trade_res, 
      gene = select_genes, 
      nPoints = 100, 
      tidy = FALSE)

  U_seq <- seq(min(field_df$U),max(field_df$U),length.out = 100) %>% as.data.frame()
  rownames(U_seq) <- colnames(yhatSmooth)
  if(is.null(pot_col)){
    col_fun <- circlize::colorRamp2(c(min(field_df$U), 0, max(field_df$U)), c("dodgerblue3","lightgray","salmon3"))
  }else(
    col_fun <- pot_col
  )
  
  ha <- HeatmapAnnotation(
      potential = anno_simple(U_seq,col = col_fun,border=TRUE,which="column"),
      show_annotation_name = FALSE
      )

  legend <- list(
      title = "Relative Expression",
      legend_height = unit(6, "cm"),
      title_position = "leftcenter-rot",
      title_gp = gpar(fontsize = 12),
      border="black"
        #legend_height = unit(20, "mm")
        )

  yhat_scaled <- yhatSmooth %>% t %>% scale() %>% t
  set.seed(42) #set.seed for reproducible cluster block order
  ht <- 
    ComplexHeatmap::Heatmap(yhat_scaled,
      row_km = n_split_row,
      show_parent_dend_line = FALSE,
      row_title = NULL,
      col=ht_col,
      cluster_rows = TRUE,cluster_columns = FALSE,
      show_row_names = show_rownames,show_column_names = FALSE,
      column_title = "Field Potential", 
      column_title_side = "bottom",
      bottom_annotation = ha,
      heatmap_legend_param = legend
    )
  main_title <- 
    ifelse(is.null(title),
      paste0(stringr::str_replace(field_df$LR_pair[1],pattern = "\\.",replacement = "->"),
            " Field Potential DE genes"),
      title
    )

  ht2 <- draw(ht, 
    column_title = main_title,
    column_title_gp = gpar(fontsize = 16))
  set.seed(NULL)
  return(ht2)
}

#************************************
#** multi layers comparison Module **
#************************************

#' aggregate multi layers C2C score
#'
#' Aggregate multi layers C2C score in one dataframe, a helper function of plot_layers_score_dot
#'
#' @param score_list merged list of C2C score result list.
#' @param label_name character vector.
#' @param merge_db merged db list.
#'
#' @return return aggregated C2C score dataframe
#' @import dplyr
#' @importFrom tidyr drop_na
#' @export
aggregate_layers_C2C <- function(score_list,label_name,merge_db){
  res <- 
    aggregate_C2C_score(
      db_C2C_score_list=score_list,
      kept_db=merge_db
    ) %>%
    tidyr::drop_na() %>% 
    dplyr::mutate(label = label_name)
  return(res)
}


#' Plot multi layers C2C score as Dot Plot
#'
#' Plot multi C2C score as Dot Plot 
#'
#' @param merge_db_C2C_list list of C2C score result list.
#' @param merge_kept_db_list list of matching kept_db.
#' @param p_val p_val threshold, default is 0.05.
#' @param LR_pair LR pair to be shown, must in format "Ligand.Receptor".
#' @param LR_family LR family to be shown.
#' @param source_use Source cluster to be shown.
#' @param target_use Target cluster to be shown.
#' @param scale logical, whether to scale the score, default is TRUE.
#'
#' @return return a dot plot
#' @importFrom purrr map2
#' @import stringr
#' @import dplyr
#' @export
plot_layers_score_dot <- function(
  merge_db_C2C_list,
  merge_kept_db_list,
  p_val=0.05,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  scale=TRUE
){
  merge_db <- 
    do.call(rbind.data.frame,merge_kept_db_list) %>%
    dplyr::distinct(.keep_all = TRUE)
  res_list <- 
    purrr::map2(.x=merge_db_C2C_list,.y=names(merge_db_C2C_list),
                .f=aggregate_layers_C2C,merge_db=merge_db)
  merge_res_df <- do.call(rbind.data.frame,res_list)
  merge_res_df$label %<>% factor(levels=names(merge_db_C2C_list))

  #** 后面就是照抄的
  temp <- merge_res_df %>% 
    prep_C2C_plot_df( 
      LR_pair=LR_pair,
      LR_family=LR_family,
      source_use=source_use,
      target_use=target_use,
      scale=scale
    )

  draw_df <- temp %>% 
    dplyr::filter(p_value<0.05) %>%
    dplyr::mutate(
        LR_pair=stringr::str_replace(.data[["id"]],pattern = "\\.",replacement = "->"),
        S2T=paste0(.data[["Source"]],"->",.data[["Target"]])
        )
  if(scale && nrow(draw_df)!=1){
    draw_df <- draw_df %>% 
      dplyr::mutate(scale_score = scale(.data[["raw_score"]],center = FALSE))
    score_title <- "scaled score"
  }else{
    if(nrow(draw_df)==1){message("Only one row in filtered data. Will not scale.")}
    draw_df <- draw_df %>% 
      dplyr::mutate(scale_score = .data[["raw_score"]])
    score_title <- "score"
  }
  plot <- 
    draw_df %>%
      ggplot(aes(x=S2T,y=LR_pair)) +
        geom_point(aes(x=S2T,y=LR_pair,color=p_value,size=scale_score))+
        scale_color_gradient(low = "firebrick3",high = "navy") +
        theme_minimal()+
        scale_x_discrete(guide = guide_axis(angle = 90))+
        labs(size=score_title) + 
        facet_wrap(~label,scales = "free_x")+
        theme(strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "white")
        )
  return(plot)
}

#** adapt for spWAVE class
#' Plot multi layers C2C score as Dot Plot for spWAVE class list
#'
#' Plot multi C2C score as Dot Plot for list of spWAVE class object
#' 
#' @param merge_db_C2C_list list of C2C score result list.
#' @param merge_kept_db_list list of matching kept_db, Default is NULL. optional.
#' @param p_val p_val threshold, default is 0.05.
#' @param LR_pair LR pair to be shown, must in format "Ligand.Receptor".
#' @param LR_family LR family to be shown.
#' @param source_use Source cluster to be shown.
#' @param target_use Target cluster to be shown.
#' @param scale logical, whether to scale the score, default is TRUE.
#' 
#' @seealso \link{plot_layers_score_dot}
#'
#' @return return a dot plot
#' @importFrom purrr map2
#' @import stringr
#' @import dplyr
#' @export
plot_layers_score_dot_spWAVE <- function(
  merge_db_C2C_list,
  merge_kept_db_list=NULL,
  p_val=0.05,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  scale=TRUE
){
  check_class <- sapply(merge_db_C2C_list,function(x) is(x,"spWAVE"))
  if(!all(check_class)){
    stop("The input contents in list must be spWAVE class!")
  }
  extract_db_field_list <- lapply(merge_db_C2C_list,function(x) x@C2C_score)
  names(extract_db_field_list) <- names(merge_db_C2C_list)
  if(is.null(merge_kept_db_list)){
    extract_kept_db_list <- lapply(merge_db_C2C_list,function(x) x@kept_db)
    names(extract_kept_db_list) <- names(merge_db_C2C_list)
  }

  plot_layers_score_dot(
    merge_db_C2C_list=extract_db_field_list,
    merge_kept_db_list=extract_kept_db_list,
    p_val=p_val,
    LR_pair=LR_pair,
    LR_family=LR_family,
    source_use=source_use,
    target_use=target_use,
    scale=scale
  )
}