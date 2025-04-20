#' Plot field projection on image
#'
#' Plot given LR field projection over image of HE，cluster,
#' expression or other quantity.
#'
#' @param arrow_df data.frame, contain barcode, coordinates, field vector, expression and other quantities.
#' Usually using the output from \code{\link{extract_LR_field_result}}.
#' @param seurat_obj Seurat object, used to extract image.
#' @param mode arrow show in each point (spot) or in grid.
#' Only support "point" or "grid", default is "grid".
#' @param point point (spot) value to show, must be in colnames of arrow_df.
#' @param point_size point (spot) size, default is 2.
#' @param show_arrow logical, whether to show arrow, default is TRUE.
#' @param arrow_sf arrow head size scale factor, smaller value means larger head size.
#' @param line_sf arrow line size scale factor, smaller value means larger line size.
#' @param grid_density grid density of arrow, smaller value means sparser arrow over image, default is 0.5.
#' @param image image type, "blank" or "HE" or "cluster", default is "blank".
#' @param arrow_color arrow color.
#' @param point_color point color of value.
#' @param arrow_normalize logical, TRUE means only show the direction without arrow length, default is FALSE.
#' @param ... Other args passing to \code{\link{SpatialPlot}} in Seurat.
#' @details we do recommend to use "grid" instead of "point".
#' @return field projection plot
#' @export
plot_field_direction <- function(
  arrow_df,
  seurat_obj,
  mode=c("grid","point"),
  point="E_strength",
  point_size=2,
  show_arrow=TRUE,
  arrow_sf=NULL,
  line_sf=NULL,
  grid_density=NULL,
  image=c("blank","HE","cluster"),
  arrow_color=NULL,
  point_color="navy",
  arrow_normalize=FALSE,
  ...
  #alpha=c(0.5,1),
  ){
  #check arg
  mode <- match.arg(mode)
  image <- match.arg(image)
  #create the data.frame for drawing

  if(mode=="grid"){
    arrow_sf <- ifelse(is.null(arrow_sf),0.2,arrow_sf)
    line_sf <- ifelse(is.null(line_sf),1,line_sf)
    grid_density <- ifelse(is.null(grid_density),0.5,grid_density)
    arrow_draw <- generate_grid_vector(arrow_df,grid_density = grid_density)
  }else{
    arrow_sf <- ifelse(is.null(arrow_sf),1,arrow_sf)
    line_sf <- ifelse(is.null(line_sf),5,line_sf)
    arrow_draw <- arrow_df
  }

  if(show_arrow){
  }else{
    arrow_sf <- 0
    line_sf <- Inf
  }

  if(is.null(arrow_color)){
    if(image=="blank"){
      arrow_color <- "firebrick3"
    }else{
      arrow_color <- "black"
    }
  }

  arrow_draw %<>%
    calc_field_strength() %>%
    optimize.arrow(scale.factor = line_sf,normalize = arrow_normalize)
  arrow_df %<>% calc_field_strength()

  if(arrow_normalize){
    line_sf <- ifelse(is.null(line_sf),50,line_sf)
    arrow_size <- line_sf/4000
  }else{
    arrow_size <- arrow_draw$E_strength/arrow_sf
  }

  #draw image in each layer
  switch(image,
    blank={
      p_arrow <-
        ggplot(data = arrow_df,aes(x=x,y=y)) +
          geom_point(
            data=arrow_df,
            aes(x=x,y=y,fill=.data[[point]]),
            size=point_size, shape = 23, stroke = 0.1
          )+
          scale_fill_gradient(low="lightgrey", high=point_color) +
          geom_segment(
            data=arrow_draw,
            aes(x=x,y=y,xend = x + Ex.u, yend = y + Ey.u),
            linewidth = 1, color = arrow_color,
            arrow = arrow(length = unit(arrow_size, "npc"))
          )+
        #scale_colour_hue(l = 45) +
          theme_classic() +
          theme(axis.text.x = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0, hjust = 1),
                axis.text.y = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0))+
          coord_equal()
    },
    HE={
      p_arrow <-
        (SpatialFeaturePlot(seurat_obj, features = NULL, alpha = c(0)) + NoLegend())+
          geom_segment(
            data=arrow_draw,
            aes(x=x,y=y,xend = x + Ex.u, yend = y + Ey.u,fill=NULL),
            linewidth = 1, color = arrow_color,
            arrow = arrow(length = unit(arrow_size, "npc"))
          )+
        #scale_colour_hue(l = 45) +
          theme_classic() +
          theme(axis.text.x = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0, hjust = 1),
                axis.text.y = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0))+
          coord_equal()
    },
    cluster={
      p_arrow <-
        (SpatialDimPlot(seurat_obj, label = FALSE,...) + NoLegend())+
          geom_segment(
            data=arrow_draw,
            aes(x=x,y=y,xend = x + Ex.u, yend = y + Ey.u,fill=NULL),
            linewidth = 1, color = arrow_color,
            arrow = arrow(length = unit(arrow_size, "npc"))
          )+
        #scale_colour_hue(l = 45) +
          theme_classic() +
          theme(axis.text.x = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0, hjust = 1),
                axis.text.y = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0))+
          coord_equal()
    }
  )
  return(p_arrow)
}

#' Another Plot field projection on image
#'
#' Plot given LR field projection over image of HE，cluster,
#' expression or other quantity, using ggquiver.
#'
#' @param arrow_df data.frame, contain barcode, coordinates, field vector, expression and other quantities.
#' Usually using the output from \code{\link{extract_LR_field_result}}.
#' @param seurat_obj Seurat object, used to extract image.
#' @param mode arrow show in each point (spot) or in grid.
#' Only support "point" or "grid", default is "grid".
#' @param point point (spot) value to show, must be in colnames of arrow_df.
#' @param point_size point (spot) size. In 10x Visium result, default is 2.
#' In centroid result, default is 0.2.
#' @param point_color point color of value. Default is navy.
#' @param hexified hexified the point display, default is FALSE.
#' @param hex_bins hexified bin number, default is NULL.
#' @param hex_binwidth hexified bin width, numeric vector giving bin width in both vertical and horizontal directions.
#' Default is c(50,50) micron. Overrides hex_bins if both set.
#' @param grid_density grid density of arrow, smaller value means sparser arrow over image, default is 0.5.
#' @param image image type, "blank" or "HE" or "cluster", default is "blank".
#' @param show_arrow logical, whether to show arrow, default is TRUE.
#' @param arrow_sf scale factor of arrow, smaller value means smaller size. Default is 2.
#' And may only effect in mode="point", because of the ggquiver behavior. See details.
#' @param arrow_linewidth the linewidth of arrow. Default: 1
#' @param arrow_color arrow color. Default: firebrick3 in blank, black in other.
#' But in 'cluster' of centroid result, the color is white to match the black background.
#' @param arrow_alpha arrow alpha. Default: 1 in blank, 0.7 in other.
#' @param arrow_shadow whether add shadow under arrow. Default: TRUE.
#' @param shadow_adj shadow position adjust. Default: 0.3.
#' @param shadow_color shadow color. Default: "#c3c3c3".
#' @param shadow_alpha shadow alpha. Default: 0.7.
#' while shadow_alpha larger than arrow_alpha,
#' shadow_alpha will be adjust by the arrow_alpha.
#' @param arrow_normalize logical, TRUE means only show the direction without arrow length, default is FALSE.
#' @param return_data default is FALSE, whether return drawing data for custom plot function.
#' @param ... Other args passing to \code{\link{SpatialPlot}} or \code{\link{ImageDimPlot}} in Seurat.
#' @details we do recommend to use "grid" instead of "point".
#' This function using ggquiver to draw the vector field projection differed from plot_field_direction().
#' ggquiver will auto resize the arrow in grid mode, but not in point mode.
# ' The `arrow_sf` is a scale factor used to fine-tuning the size of the arrows.
# ' The final size is determined by arrow_sf multiplied with a scaling magnitude,
# ' which is determined based on the magnitude of the vectors.
#'
#' @return field projection plot
#' @import ggplot2
#' @import ggquiver
#' @importFrom Seurat SpatialDimPlot ImageDimPlot SpatialFeaturePlot
#' @export
plot_field_direction2 <- function(
  arrow_df,
  seurat_obj,
  mode=c("grid","point"),
  point="E_strength",
  point_size=NULL,
  point_color="navy",
  hexified=FALSE,
  hex_bins=NULL,
  hex_binwidth=c(50,50),
  grid_density=NULL,
  image=c("blank","HE","cluster"),
  show_arrow=TRUE,
  arrow_sf=2,
  arrow_linewidth=1,
  arrow_color=NULL,
  arrow_alpha=NULL,
  arrow_shadow=TRUE,
  shadow_adj=0.3,
  shadow_color="#c3c3c3",
  shadow_alpha=0.7,
  arrow_normalize=FALSE,
  return_data=FALSE,
  ...
  #alpha=c(0.5,1),
  ){
  #check arg
  mode <- match.arg(mode)
  image <- match.arg(image)
  #create the data.frame for drawing
  #** first distinguish seurat class and use different function and default args
  if(image=="blank"){
      arrow_color <- ifelse(is.null(arrow_color),"firebrick3",arrow_color)
      arrow_alpha <- ifelse(is.null(arrow_alpha),1,arrow_alpha)
      if(is.null(point_size)){
        point_size <-
          ifelse(nrow(arrow_df)<10000,2,
                ifelse(nrow(arrow_df)<20000,1,0.2)
                )
      }
    }else{
      arrow_alpha <- ifelse(is.null(arrow_alpha),0.7,arrow_alpha)

      if(class(seurat_obj@images[[1]]) %in% c("FOV")){
        #** point_size only works in "blank", so not need to add if statement. Below is same.
        if(image == "HE"){
          warning("The centroid methods do not have HE image, using 'cluster' instead.")
          image <- "cluster"
        }
        if(image == "cluster"){
          SeuratDimPlot <- Seurat::ImageDimPlot
          arrow_color <- ifelse(is.null(arrow_color),"white",arrow_color)
        }
      }else if(class(seurat_obj@images[[1]]) %in% c("VisiumV1","VisiumV2")){
        temp <- rownames(arrow_df)
        lo_res_coord <- get_coordinates_in_plot(seurat_obj)
        arrow_df %<>% dplyr::rows_update(lo_res_coord,by='barcode')
        rownames(arrow_df) <- arrow_df$barcode
        arrow_df <- arrow_df[temp,]

        arrow_color <- ifelse(is.null(arrow_color),"black",arrow_color)
        if(image=="cluster"){
          SeuratDimPlot <- Seurat::SpatialDimPlot
        }
      }else{
        stop("Invaild Seurat Object! Please check the 'images' slot.")
      }
    }
  if(mode=="grid"){
    #arrow_sf <- ifelse(is.null(arrow_sf),0.2,arrow_sf)
    #line_sf <- ifelse(is.null(line_sf),1,line_sf)
    grid_density <- ifelse(is.null(grid_density),0.5,grid_density)
    arrow_draw <- generate_grid_vector(arrow_df,grid_density = grid_density)
  }else{
    #arrow_sf <- ifelse(is.null(arrow_sf),1,arrow_sf)
    #line_sf <- ifelse(is.null(line_sf),5,line_sf)
    arrow_draw <- arrow_df
  }

  if(show_arrow){
    arrow_sf <- autocalc_arrow_sf(arrow_draw)*arrow_sf
  }else{
    #arrow_sf <- 0
    #line_sf <- Inf
    arrow_sf=0
  }

  if(arrow_normalize){
    arrow_draw %<>%
      mutate(length=sqrt(Ex^2+Ey^2),
      Ex=Ex/length,
      Ey=Ey/length
      )
  }
  #arrow_draw %<>%
    #calc_field_strength() %>%
    #optimize.arrow(normalize = arrow_normalize)
  arrow_df %<>% calc_field_strength()

  #draw image in each layer
  #bg layer
  switch(image,
    blank={
      if(hexified){
        p_arrow <-
          ggplot(data = arrow_df,aes(x=x,y=y)) +
            stat_summary_hex(
              data=arrow_df,
              fun=mean,
              aes(x=x,y=y,z=.data[[point]]),
              bins = hex_bins,
              binwidth = hex_binwidth
            )+
            scale_fill_gradient(low="lightgrey", high=point_color) +
          #scale_colour_hue(l = 45) +
            theme_classic()
      }else if(is.numeric(arrow_df[[point]])){
        p_arrow <-
          ggplot(data = arrow_df,aes(x=x,y=y)) +
            geom_point(
              data=arrow_df,
              aes(x=x,y=y,fill=.data[[point]]),
              size=point_size, shape = 23, stroke = 0.1
            )+
            scale_fill_gradient(low="lightgrey", high=point_color) +
          #scale_colour_hue(l = 45) +
            theme_classic()
      }else{
        p_arrow <-
          ggplot(data = arrow_df,aes(x=x,y=y)) +
            geom_point(
              data=arrow_df,
              aes(x=x,y=y,colour=.data[[point]]),
              size=point_size,
            )+
            theme_classic()

        if(
          length(point_color)>=length(unique(arrow_df[[point]]))
          ){
          p_arrow <-
            p_arrow + scale_colour_manual(values=point_color)
        }else{
          warning("Insufficient number given in point_color, use default palette instead")
        }
      }
    },
    HE={
      p_arrow <-
        SpatialFeaturePlot(seurat_obj, features = NULL, alpha = c(0)) + NoLegend()
    },
    cluster={
      p_arrow <-
        SeuratDimPlot(seurat_obj,...) + NoLegend()
    }
  )

  #arrow layer
  if(show_arrow){
    if(arrow_shadow){ # add arrow shadow
      if(shadow_alpha>arrow_alpha){
        shadow_alpha <- arrow_alpha*shadow_alpha
      }
      p_arrow <-
        p_arrow +
          geom_quiver(
            data=arrow_draw,
            aes(x=x+shadow_adj,y=y+shadow_adj,
                u=arrow_sf*Ex+shadow_adj,v=arrow_sf*Ey+shadow_adj,fill=NULL),
            linewidth = arrow_linewidth, color = shadow_color,alpha=shadow_alpha
          )
    }
    p_arrow <-
      p_arrow +
        geom_quiver(
          data=arrow_draw,
          aes(x=x,y=y,u=arrow_sf*Ex, v=arrow_sf*Ey,fill=NULL),
          linewidth = arrow_linewidth, color = arrow_color,alpha=arrow_alpha
        )
  }

  #theme adjust
  p_arrow <-
    p_arrow +
      theme(axis.text.x = element_text(face = "bold", color = "black",
                                      size = 12, angle = 0, hjust = 1),
            axis.text.y = element_text(face = "bold", color = "black",
                                            size = 12, angle = 0))+
      coord_equal()
  if(return_data){
    return(
      list(p=p_arrow,
        spot=arrow_df,
        quiver=arrow_draw)
      )
  }
  else{
    return(p_arrow)
  }
}

#' Cluster levels interaction scores dot plot
#'
#' Dot Plot to display cluster levels interaction scores from database calc result
#'
#' @param db_C2C_score_list C2C score result list.
#' @param kept_db The database used in results list.
#' @param p_val p value cut off, default is 0.05.
#' @param LR_pair LR pair selected to display, must be the format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param LR_family LR family selected to display.
#' @param source_use Ligand source clusters selected to display.
#' @param target_use Receptor target clusters selected to display.
#' @param scale logical, whether to scale the score, default is TRUE.
#' @import dplyr
#' @import stringr
#' @return return a ggplot2 object plot
#' @export
plot_db_score_dot_impl <- function(
  db_C2C_score_list,
  kept_db,
  p_val=0.05,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  scale=TRUE
){
  aggregate_df <-
    aggregate_C2C_score(
      db_C2C_score_list=db_C2C_score_list,
      kept_db=kept_db
    ) %>%
    prep_C2C_plot_df(
        LR_pair=LR_pair,
        LR_family=LR_family,
        source_use=source_use,
        target_use=target_use,
        scale=scale
    )

  draw_df <-
    aggregate_df %>%
    dplyr::filter(p_value<p_val) %>%
      dplyr::mutate(
        LR_pair=stringr::str_replace(.data[["id"]],pattern = "\\.",replacement = "->"),
        S2T=paste0(.data[["Source"]],"->",.data[["Target"]])
        )
  if(scale && nrow(draw_df)!=1){
    draw_df <- draw_df %>%
      dplyr::mutate(scale_score = scale(.data[["raw_score"]]))
    score_title <- "scaled score"
  }else{
    if(nrow(draw_df)==1){message("Only one row in filtered data. Will not scale.")}
    draw_df <- draw_df %>%
      dplyr::mutate(scale_score = .data[["raw_score"]])
    score_title <- "score"
  }

  p_dot <-
  draw_df %>%
    ggplot(aes(x=S2T,y=LR_pair)) +
      geom_point(aes(x=S2T,y=LR_pair,color=p_value,size=scale_score))+
      scale_color_gradient(low = "firebrick3",high = "navy") +
      theme_minimal()+
      scale_x_discrete(guide = guide_axis(angle = 90))+
      labs(size=score_title)
  return(p_dot)
}

#' Material Design 2 color palette picker
#'
#' Material Design 2 color palette picker from ggsci
#'
#' @param n number of colors, maximum is 37.
#'
#' @return return color palette
#' @import ggsci
#' @export
MD2_color_picker <- function(n){
  color_order <- c("red", "pink", "purple", "deep-purple", "indigo", "blue", "light-blue",
    "cyan", "teal", "green", "light-green", "lime", "yellow", "amber", "orange",
    "deep-orange", "brown","blue-grey")
  if(n<=19){
    #color_order <- c(color_order,"grey")
    MD_mix_pal <-
      sapply(color_order,function(color){
      p_color <- pal_material(color)(10)[6]
      return(p_color)
    })
    c(MD_mix_pal,"9E9E9E")[1:n]
  }else if(n<=37){
    pal1 <-
      sapply(color_order,function(color){
      p_color <- pal_material(color)(10)[6]
      return(p_color)
    })
    color_order <- c(color_order,"grey")
    pal2 <-
      sapply(color_order,function(color){
      p_color <- pal_material(color)(10)[9]
      return(p_color)
    })
    c(pal1,pal2)[1:n]
  }else{
    stop("Too many colors!")
  }
}



#' Core function of chord plot series functions
#'
#' Core function of chord plot series functions
#'
#' @param pic_df dataframe contain interactiion with source , target, and score.
#' @param grid_col named vector, the values are colors and the names are cluster.
#' @param chord_order cluster order.
#' @param chord_group cluster names.
#' @param main_title plot title.
#' @importFrom graphics strwidth
#' @importFrom circlize circos.clear circos.par chordDiagram circos.text circos.track get.cell.meta.data highlight.sector
#' @import stringr
#' @export
chord_cluster_core_function <- function(
  pic_df,grid_col,chord_order,chord_group,main_title
){
  circos.clear()
  circos.par(start.degree = 175)
  chordDiagram(pic_df[,c(1:3)],
      grid.col = grid_col,
      #col=link_color,
      directional = 1,big.gap = 10,
      direction.type = c("diffHeight", "arrows"),
      group = chord_group,
      order=chord_order,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = max(graphics::strwidth(chord_order))),
      link.arr.type = "big.arrow",
      scale = FALSE
      )

  circos.track(track.index = 2, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      xplot = get.cell.meta.data("xplot")
      ylim = get.cell.meta.data("ylim")
      sector.name = stringr::str_split(get.cell.meta.data("sector.index"),pattern = "@",simplify = TRUE)[,2]
      circos.text(mean(xlim), ylim[2], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5),cex = 1)
    }, bg.border = NA) # here set bg.border to NA is important
  title(main_title)

  highlight.sector(chord_order[stringr::str_detect(chord_order,pattern = "^S")], track.index = 1,col = NA,
      text = "Sender",cex = 1.2, facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  highlight.sector(chord_order[stringr::str_detect(chord_order,pattern = "^R")], track.index = 1,col = NA,
      text = "Receiver",cex = 1.2,facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  circos.clear()
}



#' Chordigram plot for cluster level Ligand-receptor pairs
#'
#' Chordigram plot for cluster level Ligand-receptor pairs
#'
#' @param db_C2C_score_list C2C score result list.
#' @param kept_db The database used in results list.
#' @param cluster_color named vector, the values are colors and the names are cluster.
#' @param LR_pair LR pair selected to display, must be the format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param LR_family LR family selected to display.
#' @param source_use Ligand source clusters selected to display.
#' @param target_use Receptor target clusters selected to display.
#' @param title character, plot title.
#' Default is NULL and can be generated automatically when single LR pair or family given.
#' @param scale logical, whether to scale the arc of each chord in plot, default is FALSE.
#'
#' @return plot using based on circlize package
#'
#' @import stringr
#' @import dplyr
#' @export
plot_LR_cluster_chord_impl <- function(
  db_C2C_score_list,
  kept_db,
  #C2C_score_df,
  cluster_color=NULL,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  title=NULL,
  scale=FALSE
){
  #** assign color to all clusters
  C2C_score_df <-
    aggregate_C2C_score(db_C2C_score_list,kept_db)

  #** get the color of each cluster
  if(is.null(cluster_color)){
    grid_col <- assign_cluster_color(C2C_score_df)

  }else{
    message("The clusters in plot will follow the order of cluster_color names.")
    clu_color <- check_cluster_color(cluster_color,C2C_score_df,source_use,target_use)
    rev_color <- rev(clu_color)
    names(rev_color) <- paste0(rep("R@",length(clu_color)),rev(names(clu_color)))
    names(clu_color) <- paste0(rep("S@",length(clu_color)),names(clu_color))
    grid_col <- c(clu_color,rev_color)

  }
  chord_order <- names(grid_col)

  #** filter the data to draw plot
  pic_df <-
    C2C_score_df %>%
      prep_C2C_plot_df(
          LR_pair=LR_pair,
          LR_family=LR_family,
          source_use=source_use,
          target_use=target_use,
          scale=scale
      ) %>%
      dplyr::filter(p_value<0.05) %>%
      #dplyr::filter(Source!=Target) %>%
      dplyr::mutate(Source = paste0("S@",Source),Target = paste0("R@",Target))

  temp <- which(chord_order %in% c(pic_df$Source,pic_df$Target))
  chord_order <- chord_order[temp]
  grid_col <- grid_col[temp]
  chord_group <- stringr::str_split(chord_order,pattern = "@",simplify = TRUE)[,1]
  names(chord_group) <- chord_order

  title0 <- NULL
  if(!is.null(LR_family)){title0 <- paste0(LR_family," family signaling")}
  if(length(LR_pair)==1){
    title0 <- stringr::str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }
  if(is.null(title)){title <- title0}

    chord_cluster_core_function(
      pic_df,
      grid_col,
      chord_order,chord_group,title
  )
}



#' Core function of chord plot for genes
#'
#' Core function of chord plot for genes
#'
#' @param pic_df dataframe contain interactiion with source , target, and score.
#' @param grid_col named vector, color for arc, the values are colors and the names are cluster, usually same with legend_col
#' @param legend_col named vector, color for legend, the values are colors and the names are cluster, usually same with grid_col.
#' @param chord_order cluster order.
#' @param chord_group cluster names.
#' @param main_title plot title.
#' @import ComplexHeatmap
#' @importFrom circlize circos.clear circos.par chordDiagram circos.text circos.track get.cell.meta.data highlight.sector
#' @import stringr
#' @importFrom graphics title strwidth
#' @export
chord_gene_core_function <- function(
  pic_df,grid_col,legend_col,chord_order,chord_group,main_title
){
  circos.clear()
  circos.par(start.degree = 176)
  chordDiagram(pic_df[,c("lig","rec","raw_score")],
      grid.col = grid_col,
      #col=link_color,
      directional = 1,big.gap = 10,
      direction.type = c("diffHeight", "arrows"),
      group = chord_group,
      order=chord_order,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = max(graphics::strwidth(chord_order))),
      link.arr.type = "big.arrow",
      scale = FALSE
      )

  circos.track(track.index = 2, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      xplot = get.cell.meta.data("xplot")
      ylim = get.cell.meta.data("ylim")
      temp = stringr::str_split_fixed(get.cell.meta.data("sector.index"),pattern = "[@\\.]",n=2)
      sector.name = temp[,2]
      circos.text(mean(xlim), ylim[2], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5),cex = 1)
    }, bg.border = NA) # here set bg.border to NA is important

  highlight.sector(chord_order[stringr::str_detect(chord_order,pattern = "^S")], track.index = 1,col = NA,
      text = "Ligands",cex = 1.2, facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  highlight.sector(chord_order[stringr::str_detect(chord_order,pattern = "^R")], track.index = 1,col = NA,
      text = "Receptors",cex = 1.2,facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)

  #legend <- ComplexHeatmap::Legend(at = names(legend_col), type = "grid", legend_gp = grid::gpar(fill = legend_col), title = "Cluster")
  #ComplexHeatmap::draw(legend, x = unit(1, "npc")-unit(20, "mm"), just = c("right"))

  graphics::title(main_title)

  circos.clear()
}



#' Chordigram plot for genes in Ligand-receptor pairs
#'
#' Chordigram plot for genes in Ligand-receptor pairs
#'
#' @param db_C2C_score_list C2C score result list.
#' @param kept_db The database used in results list.
#' @param LR_color named vector, the values are colors and the names are ligand and receptor.
#' @param LR_pair LR pair selected to display, must be the format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param LR_family LR family selected to display.
#' @param source_use Ligand source clusters selected to display.
#' @param target_use Receptor target clusters selected to display.
#' @param title character, plot title.
#' Default is NULL and can be generated automatically when single LR pair or family given.
#' @param scale logical, whether to scale the arc of each chord in plot, default is FALSE.
#'
#' @return plot using based on circlize package
#' @import dplyr
#' @import stringr
#' @export
plot_LR_gene_chord_impl <- function(
  db_C2C_score_list,
  kept_db,
  #C2C_score_df,
  LR_color=NULL,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  title=NULL,
  scale=FALSE
){
  C2C_score_df <-
    aggregate_C2C_score(db_C2C_score_list,kept_db)
  #** 这里需要添加一个转换，把cluster名字里面的.转换成下划线_,不然会干扰到后面LR判断
  message("The '.' in the 'cluster' will be transformed to '_' in case of conflicting with LR")
  #** filter the data to draw plot
  pic_df <-
    C2C_score_df %>%
      prep_C2C_plot_df(
          LR_pair=LR_pair,
          LR_family=LR_family,
          source_use=source_use,
          target_use=target_use,
          scale=scale
      ) %>%
      dplyr::filter(p_value<0.05) %>%
      #dplyr::filter(Source!=Target) %>%
      #dplyr::mutate(Source = paste0("S@",Source),Target = paste0("R@",Target)) %>%
      #dplyr::mutate(Source=stringr::str_replace_all(Source,pattern = "\\.",replacement = "_"),
      #      Target=stringr::str_replace_all(Target,pattern = "\\.",replacement = "_")) %>%
      #dplyr::mutate(Source = paste0(Source,"_",id),Target = paste0(Target,"_",id))
      dplyr::mutate("lig"=stringr::str_split(id,pattern = "\\.",simplify = TRUE)[,1],
            "rec"=stringr::str_split(id,pattern = "\\.",simplify = TRUE)[,2]) %>%
      dplyr::mutate(lig=paste0("S@",lig),rec=paste0("R@",rec)) #%>%
      # dplyr::mutate(source_lig=paste0(Source,".",lig),
      #       target_rec=paste0(Target,".",rec))


  #** get the color of each LR
  if(is.null(LR_color)){
    LR_color <- MD2_color_picker(length(c(pic_df$lig,pic_df$rec) %>% unique()))
    names(LR_color) <- c(pic_df$lig,pic_df$rec) %>% unique()
  }else{
    #message("The clusters in plot will follow the order of cluster_color names.")
    temp1 <- which(paste0("S@",names(LR_color)) %in% pic_df$lig)
    temp2 <- which(paste0("R@",names(LR_color)) %in% pic_df$rec)
    names(LR_color)[temp1] <- paste0("S@",names(LR_color)[temp1])
    names(LR_color)[temp2] <- paste0("R@",names(LR_color)[temp2])
    LR_color <- c(LR_color[temp1],LR_color[temp2])
  }


  #** assign the color of each cluster and ligand/receptor
  # color_df <- clu_color %>% as.data.frame()
  # colnames(color_df) <- "color"
  # color_df$clu <- rownames(color_df)
  # LR_clu <- c(pic_df$source_lig,pic_df$target_rec) %>% unique() %>%
  #   stringr::str_split_fixed(pattern="[@\\.]",n=3) %>% as.data.frame()
  # colnames(LR_clu) <- c("SR","clu","LR")

  # LR_clu %<>% dplyr::left_join(color_df,by="clu")
  # grid_col <- LR_clu$color
  chord_order <-
    #names(grid_col) <-
      c(pic_df$lig,pic_df$rec) %>% unique()


  #** filter out not used color,and this arg is used in legend
  #clu_color_use <- clu_color[LR_clu$clu %>% unique()]

  #** filter out not used color,and those args is used in chord plot
  #temp <- which(chord_order %in% c(pic_df$lig,pic_df$rec))
  #chord_order <- chord_order[temp]
  grid_col <- LR_color[chord_order]
  chord_group <- stringr::str_split(chord_order,pattern = "@",simplify = TRUE)[,1]
  names(chord_group) <- chord_order

  #** set title
  title0 <- NULL
  if(!is.null(LR_family)){title0 <- paste0(LR_family," family signaling")}
  if(length(LR_pair)==1){
    title0 <- stringr::str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }
  if(is.null(title)){title <- title0}

  chord_gene_core_function(
      pic_df,
      grid_col,legend_col=NULL,
      chord_order,chord_group,title
  )
}

#' Heatmap of ligand-receptor interaction
#'
#' Heatmap of ligand-receptor interaction in clusters level
#'
#' @param db_C2C_score_list C2C score result list.
#' @param kept_db The database used in results list.
#' @param cluster_color named vector, the values are colors and the names are cluster.
#' @param LR_pair LR pair selected to display, must be the format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param LR_family LR family selected to display.
#' @param source_use Ligand source clusters selected to display.
#' @param target_use Receptor target clusters selected to display.
#' @param title character, plot title.
#' Default is NULL and can be generated automatically when single LR pair or family given.
#' @param method_use Statistical method to display in heatmap. Deault is "count". Only support "count" and "strength",
#' refelecting the number of interactions or the sum of strength in C2C scores, respectively.
#' @param normalize logical, whether to normalize the score in heatmap, default is FALSE.
#' @param ht_col heatmap tiles color. Support sequential palettes in \code{\link{RColorBrewer}} or single color.
#' Default is "Reds" in the palettes.
#'
#' @return heatmap using based on ComplexHeatmap package.
#' @import ComplexHeatmap
#' @import grDevices
#' @import RColorBrewer
#' @import grid
#' @import stringr
#' @export
plot_LR_cluster_heatmap_impl <- function(
  db_C2C_score_list,
  kept_db,
  #C2C_score_df,
  cluster_color=NULL,
  LR_pair=NULL,
  LR_family=NULL,
  source_use=NULL,
  target_use=NULL,
  title=NULL,
  method_use=c("count","strength"),
  normalize=FALSE,
  ht_col=NULL
){

  C2C_score_df <-
      aggregate_C2C_score(db_C2C_score_list,kept_db)
  if(is.null(cluster_color)){
    clu_col <- assign_clu_col_lite(C2C_score_df)
  }else{
    clu_col <- check_cluster_color(cluster_color,C2C_score_df,source_use,target_use)
  }

  method_use <- match.arg(method_use)

  ht_mat <-
    prep_plot_matrix(
      C2C_score_df,
      LR_pair=LR_pair,
      LR_family=LR_family,
      source_use=source_use,
      target_use=target_use,
      method_use=method_use
    )

  if(normalize){
    ht_mat <- log2(ht_mat+1)
  }

  #** set heatmap color
  ht_col <- ifelse(is.null(ht_col),"Reds",ht_col)

  avail_color <-
    c("Blues","BuGn","BuPu","GnBu","Greens","Greys","Oranges",
      "OrRd","PuBu","PuBuGn","PuRd","Purples","RdPu","Reds",
      "YlGn","YlGnBu","YlOrBr","YlOrRd")
  if(ht_col %in% avail_color){
    ht_col_use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = ht_col)))(100)
  }else if(is.color(ht_col)){
    ht_col_use = grDevices::colorRampPalette(c("white", ht_col))(100)
  }else{
    warning("Do not detect valid color, use default color 'Reds' instead.")
    ht_col_use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = "Reds")))(100)
  }

  #** decoration of heatmap
  df_col <- data.frame(group = colnames(ht_mat))
  rownames(df_col) <- colnames(ht_mat)
  df_row <- data.frame(group = rownames(ht_mat))
  rownames(df_row) <- rownames(ht_mat)

  col_col <- clu_col[rownames(df_col)]
  row_col <- clu_col[rownames(df_row)]

  col_annotation <-
    HeatmapAnnotation(df = df_col, col = list(group = col_col),
      which = "column",
      show_legend = FALSE, show_annotation_name = FALSE,
      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <-
    HeatmapAnnotation(df = df_row, col = list(group = row_col),
      which = "row",
      show_legend = FALSE, show_annotation_name = FALSE,
      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 <- rowAnnotation(
    Strength = anno_barplot(rowSums(abs(ht_mat),na.rm=TRUE),
    border = FALSE,
    gp = gpar(fill = row_col,col=row_col)),
    show_annotation_name = FALSE)

  ha2 <- HeatmapAnnotation(
      Strength = anno_barplot(colSums(abs(ht_mat),na.rm=TRUE),
      border = FALSE,
      gp = gpar(fill = col_col, col=col_col)),
      show_annotation_name = FALSE)

  legend <- list(
      title = 
        paste0(ifelse(normalize,"Scale ",""),
          ifelse(method_use=="count","count","strength")
          )
        #legend_height = unit(20, "mm")
        )

  ht <- ComplexHeatmap::Heatmap(ht_mat,
    col= ht_col_use,
    na_col = "white",
    cluster_rows = FALSE,cluster_columns = FALSE,
    bottom_annotation = col_annotation, left_annotation =row_annotation,
    top_annotation = ha2, right_annotation = ha1,
    column_title = "Target",column_names_rot = 90,column_title_side = "bottom",
    row_title = "Sources",row_title_rot = 90,row_names_side = "left",
    heatmap_legend_param = legend
  )

  #** set title
  suffix_title <- NULL
  if(!is.null(LR_family)){suffix_title <- paste0(LR_family," Family Signaling")}
  if(length(LR_pair)==1){
    suffix_title <- stringr::str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }

  main_title <-
    ifelse(is.null(title),
      paste0(ifelse(method_use=="count","Number","Strength"),
            " of ",
            ifelse(is.null(suffix_title),"Interactions",suffix_title)),
      title
    )

  ht2 <- draw(ht,
    column_title = main_title,
    column_title_gp = gpar(fontsize = 16))
  #return(ht2)
}


#** adopted from cellchat netVisual_circle

#' FUNCTION_TITLE
#'
#' FUNCTION_DESCRIPTION
#'
#' @param igraph_g DESCRIPTION.
#' @param color.use Colors represent different cell groups
#' @param main_title character, plot title
#' @param weight.scale whether scale the weight
#' @param vertex.weight The weight of vertex: either a scale value or a vector
#' @param vertex.weight.max the maximum weight of vertex; defualt = max(vertex.weight)
#' @param vertex.size.max the maximum vertex size for visualization
#' @param vertex.label.cex The label size of vertex
#' @param vertex.label.color The color of label for vertex
#' @param edge.weight.max the maximum weight of edge; defualt = max(net)
#' @param edge.width.max The maximum edge width for visualization
#' @param label.edge Whether or not shows the label of edges
#' @param alpha.edge the transprency of edge
#' @param edge.label.color The color for single arrow
#' @param edge.label.cex The size of label for arrows
#' @param edge.curved Specifies whether to draw curved edges, or not.
#' @param shape The shape of the vertex
#' @param layout The layout specification
#' @param margin The amount of empty space around the plot
#' @param vertex.size The size of vertex
#' @param arrow.width The width of arrows
#' @param arrow.size the size of arrow
#' @param text.x,text.y the x- and y-coordinates to add the text
#'
#' @importFrom igraph ends E V layout_ in_circle
#' @importFrom scales rescale
#' @details adopted from cellchat netVisual_circle function. More detail please refer to CellChat netVisual_circle
#' @return return net plot based on igraph
#' @export
plot_net_core_function <- function(igraph_g,color.use,main_title=NULL,
  weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
  edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
  edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.5, vertex.size = NULL,
  arrow.width=1,arrow.size = 0.2,
  text.x = 0, text.y = 1.5
){
  g <- igraph_g
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coord<-igraph::layout_(g,layout)
  if(nrow(coord)!=1){
    coord_scale=scale(coord)
  }else{
    coord_scale<-coord
  }

  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }

  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

  loop.angle<-ifelse(coord_scale[igraph::V(g),1]>0,-atan(coord_scale[igraph::V(g),2]/coord_scale[igraph::V(g),1]),pi-atan(coord_scale[igraph::V(g),2]/coord_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<- 0.3+edge.width.max*igraph::E(g)$weight
  }

  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-(igraph::E(g)$width/10)-0.3
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))

  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,main=main_title,
        edge.curved=edge.curved,vertex.shape=shape,layout=coord_scale,margin=margin, vertex.label.dist=label.dist,
        vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
}


#' Plot clusters level interactions in network
#'
#' Plot clusters level interactions in network with number or interaction stregth.
#'
#' @param db_C2C_score_list C2C score result list.
#' @param kept_db database used in results list.
#' @param cluster_color named vector, the values are colors and the names are cluster.
#' @param LR_pair LR pair selected to display, must be the format of Ligand.Receptor, e.g. "FGF1.FGFR1", "TGFB1.TGFBR1_TGFBR2".
#' @param LR_family LR family selected to display.
#' @param cluster_use cluster selected to display.
#' @param title character, plot title.
#' @param method_use Statistical method to display in heatmap. Deault is "count". Only support "count" and "strength",
#' refelecting the number of interactions or the sum of strength in C2C scores, respectively.
#' @param mat_scale logical, whether to scale score matrix, default is FALSE.
#' Differ from weight.scale, see details.
#' @param weight.scale logical, whether to scale edge weight and refelecting in plot, default is FALSE.
#' Differ from mat_scale, see details.
#' @param ... args passing to \code{\link{plot_net_core_function}} which adopted from CellChat netVisual_circle function.
#' See details.
#'
#' @details This funcion is wrapper and adopted from CellChat function netVisual_circle.
#' mat_scale and weight.scale are different argments,
#' The mat_scale is used to scale score matrix, and it also refelects in plot but may not be proper visualization.
#' The weight.scale is used to scale edge weight consdering the edge maximum length in plot.
#' Commonly, set mat_scale = TRUE and weight.scale = FALSE is enough to get proper visualization.
#' However, if the lines in the plot are still too large or too small, we recommond to set both mat_scale and weight.scale = TRUE.
#' @return network plot based on igraph package.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats sd
#' @import stringr
#' @export
plot_LR_cluster_net_impl <- function(
  db_C2C_score_list,
  kept_db,
  cluster_color=NULL,
  #C2C_score_df,
  LR_pair=NULL,
  LR_family=NULL,
  cluster_use=NULL,
  title=NULL,
  method_use=c("count","strength"),
  mat_scale=FALSE,
  weight.scale=FALSE,
  ...
){
  C2C_score_df <-
      aggregate_C2C_score(db_C2C_score_list,kept_db)
  if(is.null(cluster_color)){
    clu_col <- assign_clu_col_lite(C2C_score_df)
  }else{
    clu_col <- check_cluster_color(cluster_color,C2C_score_df,source_use=cluster_use,target_use=cluster_use)
  }

  method_use <- match.arg(method_use)

  ht_mat <-
    prep_plot_matrix(
      C2C_score_df,
      LR_pair=LR_pair,
      LR_family=LR_family,
      source_use=cluster_use,
      target_use=cluster_use,
      method_use=method_use
    )
  clu_col <- clu_col[rownames(ht_mat)]

  #** scale matrix
  if(mat_scale){
    #ht_mat_mean <- mean(ht_mat)
    ht_mat_std <- stats::sd(ht_mat)
    nor_mat <- ht_mat/ht_mat_std
    ht_mat <- log2(nor_mat+1)
  }

  #** set title
  suffix_title <- NULL
  if(!is.null(LR_family)){suffix_title <- paste0(LR_family," Family Signaling")}
  if(length(LR_pair)==1){
    suffix_title <- stringr::str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }

  main_title <-
    ifelse(is.null(title),
      paste0(ifelse(method_use=="count","Number","Strength"),
            " of ",
            ifelse(is.null(suffix_title),"Interactions",suffix_title)),
      title
    )
  g <- igraph::graph_from_adjacency_matrix(ht_mat, mode = "directed", weighted = T)
  plot_net_core_function(g,
    color.use=clu_col,
    main_title = main_title,
    weight.scale = weight.scale,
    ...)
}
