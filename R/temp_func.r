
calc_component_along_field <- function(source_field, direct_field,method=c("comp","cos")){
  Ta <- source_field$Ex*direct_field$Ex
  Tb <- source_field$Ey*direct_field$Ey
  temp <- (Ta+Tb)
  method <- match.arg(method)
  temp <- temp/direct_field$Enorm
  if(method=="comp"){
    return(temp)
  }else{
    source_field$Enorm <- sqrt(source_field$Ex^2+source_field$Ey^2)
    return(temp/source_field$Enorm)
  }
}

extract_LR_family_member <- function(kept_db,family){
  kept_db %>% 
    dplyr::filter(Family==family) %>%
    dplyr::pull(id)
}


calc_family_scaler_component <- function(db_field_result,kept_db,LR_family){
  LR_pair_field_list <- db_field_result$LR_pair_field
  LR_family_field <- db_field_result$LR_family_field_list[[LR_family]]

  LR_family_field$Enorm <- sqrt(LR_family_field$Ex^2+LR_family_field$Ey^2)

  family_member <- extract_LR_family_member(kept_db,LR_family)
  family_scaler_comp <- list()
  family_scaler_comp <- 
    lapply(family_member,function(LR_pair){
      calc_component_along_field(
        LR_pair_field_list[[LR_pair]],
        LR_family_field)
    })
  temp_df <- Reduce("cbind",family_scaler_comp)
  colnames(temp_df) <- family_member
  temp_df <- cbind.data.frame(LR_family_field,temp_df)
  return(temp_df)
}

calc_LR_set_field <- function(
  single_mole_field_list,
  L_genes,R_genes
){
  LR_field <- 
    calc_LR_pair_vec(single_mole_field_list,L_genes,R_genes)
  LR_expr_list <-  single_mole_field_list[c(L_genes,R_genes)]
  LR_expr_df <- do.call(cbind,lapply(LR_expr_list,function(x) x[,4,F]))

  merge_df <- cbind.data.frame(LR_field,LR_expr_df)
  merge_df %<>% 
    mutate(Rel_LR_Exp=
      rowSums(across(all_of(L_genes)))-
      rowSums(across(all_of(R_genes)))
    )
  return(merge_df)
}

calc_set_scaler_component <- function(
  LR_set_field,LR_field_list
){

  LR_set_field$Enorm <- sqrt(LR_set_field$Ex^2+LR_set_field$Ey^2)
  if(is.null(names(LR_field_list))){
    set_member <- sapply(LR_field_list,function(x){
      x$LR_pair[1]
    })
    names(LR_field_list) <- set_member
  }else{
    set_member <- names(LR_field_list)
  }

  set_scaler_comp <- list()
  set_scaler_comp <- 
    lapply(set_member,function(LR_pair){
      calc_component_along_field(
        LR_field_list[[LR_pair]],
        LR_set_field)
    })
  temp_df <- Reduce("cbind",set_scaler_comp)
  colnames(temp_df) <- paste0(set_member,"_comp")
  merge_df <- cbind.data.frame(LR_set_field,temp_df)
  
  set_scaler_comp <- 
    lapply(set_member,function(LR_pair){
      calc_component_along_field(
        LR_field_list[[LR_pair]],
        LR_set_field,
        method="cos")
    })
  temp_df <- Reduce("cbind",set_scaler_comp)
  colnames(temp_df) <- paste0(set_member,"_cos")
  merge_df <- cbind.data.frame(merge_df,temp_df)

  return(merge_df)
}

calc_set_scaler_component <- function(
  
){

}


merge_spwave_core <- function(spwave1,spwave2){
  merged_obj <- methods::new(Class="spWAVE")
  merged_obj@expr_raw <- rbind(spwave1@expr_raw,spwave2@expr_raw)
  merged_obj@expr_complex <- cbind(spwave1@expr_complex,spwave2@expr_complex)
  merged_obj@kept_db <- rbind(spwave1@kept_db,spwave2@kept_db)

  merged_obj@coord <- spwave1@coord
  merged_obj@cluster_info <- spwave1@cluster_info

  merged_obj@meta_coord <- spwave1@meta_coord
  merged_obj@meta_coord_clu <- spwave1@meta_coord_clu
  merged_obj@meta_complex <- cbind(spwave1@meta_complex,spwave2@meta_complex)

  merged_obj@single_mole_field <-
    c(spwave1@single_mole_field,spwave2@single_mole_field)
  merged_obj@LR_pair_field <-
    c(spwave1@LR_pair_field,spwave2@LR_pair_field)
  merged_obj@LR_family_field <-
    c(spwave1@LR_family_field,spwave2@LR_family_field)

  return(merged_obj)
}

update_spWAVE_object <- function(obj){
  message("Update single mole field")
  obj@single_mole_field %<>% lapply(function(x){
    x <- x[,4:7] %>% as.matrix() %>% as("dgCMatrix")
    return(x)
  })
  message("Update LR pair field")
  obj@LR_pair_field %<>% lapply(function(x){
    x <- x[,4:6] %>% as.matrix() %>% as("dgCMatrix")
    return(x)
  })
  message("Update LR Family field")
  obj@LR_family_field %<>% lapply(function(x){
    x <- x[,4:6] %>% as.matrix() %>% as("dgCMatrix")
    return(x)
  })
  return(obj)
}

theme_map_grid <- function(){
  theme_minimal() %+replace%
    theme(
      plot.background = element_rect(fill = "#FAFAFA", color = NA),
      legend.background = element_rect(fill = "#FAFAFA", color = NA),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      complete = TRUE
    )
}

coord_map_grid <- function(
  x_step = 100, y_step = 100, 
  x_expand = c(0, 50), y_expand = c(0, 0)){
    list(
      coord_fixed(),
      scale_x_continuous(
        breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = x_step),
        expand = x_expand
      ),
      scale_y_continuous(
        breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = y_step),
        expand = y_expand
      )
    )
}
