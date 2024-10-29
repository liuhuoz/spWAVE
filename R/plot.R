#' Plot field projection on image
#'
#' Plot given LR field projection over image of HE，cluster,
#' expression or other quantity.
#'
#' @param arrow_df data.frame, contain barcode, coordinates, field vector, expression and other quantities. 
#' Usually using the output from \code{\link{extract_LR_field_result}}.
#' @param seurat_obj Seurat object, used to extract image.
#' @param mode arrow show in each point (spot) or in grid. 
#' Only support "point" or "grid", default is "point".
#' @param point point (spot) value to show, must be in colnames of arrow_df.
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
#' @examples
#' # plot field projection only using arrow_df without seurat object.
#' # And point fill with expression of Sox2
#' plot_field_direction(arrow_df=brain_smad3_str, point = "Sox2", mode='grid')
#' 
#' # plot field projection over HE image
#' # HE and cluster will be extract from seurat object
#' plot_field_direction(mouse_fgfr_str,seurat_obj = brain, image="HE",mode='grid')
#' 
#' # plot field projection over HE image with passing args to SpatialPlot of Seurat
#' plot_field_direction(hgin_IGF,lze22_hgin,
#'   pt.size.factor=100,alpha=0.7, #pt.size.factor and alpha are passed to SpatialPlot
#'   image="cluster",mode="grid",
#'   show_arrow = T,arrow_color = "cyan2")
#'  
#' # plot field projection over cluster image with custom arrow color
#' # cluster info was extracted from seurat_obj@active.ident
#' plot_field_direction(mouse_fgfr_str,seurat_obj = brain,
#'   image="cluster",mode='grid',arrow_color="cyan2")
#' 
#' # plot field with custom point fill color indicating quantity in spot
#' plot_field_direction(nor_IGF,lze22_nor,point = "E_strength",mode="grid",
#'   line_sf = 0.3,arrow_sf = 0.05,
#'   point_color="salmon3",arrow_color = "black")
#' 
#' # plot field with custom triple colors indicating quantity in spot
#' plot_field_direction(temp,lze22_nor,point = "Rel_LR_Exp",mode="grid",
#'   line_sf = 0.3,arrow_sf = 0.05,arrow_color = "firebrick4",show_arrow = F) +
#' scale_fill_gradient2(low="dodgerblue3", mid="lightgray",high="salmon3")
#' 
plot_field_direction <- function(
  arrow_df,
  seurat_obj,
  mode=c("point","grid"),
  point="E_strength",
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
            size=3, shape = 23, stroke = 0.1
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
#' Only support "point" or "grid", default is "point".
#' @param point point (spot) value to show, must be in colnames of arrow_df.
#' @param show_arrow logical, whether to show arrow, default is TRUE.
#' @param scale_factor scale factor of arrow, smaller value means smaller size. Default is 100.
#' And may only effect in mode="point", because of the ggquiver behavior. See details. 
#' @param grid_density grid density of arrow, smaller value means sparser arrow over image, default is 0.5.
#' @param image image type, "blank" or "HE" or "cluster", default is "blank".
#' @param arrow_color arrow color.
#' @param point_color point color of value.
#' @param arrow_normalize logical, TRUE means only show the direction without arrow length, default is FALSE.
#' @param ... Other args passing to \code{\link{SpatialPlot}} in Seurat.
#' @details we do recommend to use "grid" instead of "point". 
#' This function using ggquiver to draw the vector field projection differed from plot_field_direction().
#' ggquiver will auto resize the arrow in grid mode, but not in point mode. 
#' The scale_factor is to control the arrow size in point mode, 
#' and proper value is ranged in 50~200 depend on interaction strength. 
#' @return field projection plot
#' @examples
#' # plot field projection only using arrow_df without seurat object.
#' # And point fill with expression of Sox2
#' plot_field_direction2(arrow_df=brain_smad3_str, point = "Sox2", mode='grid')
#' 
#' # plot field projection over HE image
#' # HE and cluster will be extract from seurat object
#' plot_field_direction2(mouse_fgfr_str,seurat_obj = brain, image="HE",mode='grid')
#' 
#' # plot field projection over HE image with passing args to SpatialPlot of Seurat
#' plot_field_direction2(hgin_IGF,lze22_hgin,
#'   pt.size.factor=100,alpha=0.7, #pt.size.factor and alpha are passed to SpatialPlot
#'   image="cluster",mode="grid",
#'   show_arrow = T,arrow_color = "cyan2")
#'  
#' # plot field projection over cluster image with custom arrow color
#' # cluster info was extracted from seurat_obj@active.ident
#' plot_field_direction2(mouse_fgfr_str,seurat_obj = brain,
#'   image="cluster",mode='grid',arrow_color="cyan2")
#' 
#' # plot field with custom point fill color indicating quantity in spot
#' plot_field_direction2(nor_IGF,lze22_nor,point = "E_strength",mode="grid",
#'   point_color="salmon3",arrow_color = "black")
#' 
#' # plot field with custom triple colors indicating quantity in spot
#' plot_field_direction(temp,lze22_nor,point = "Rel_LR_Exp",mode="grid",
#'   arrow_color = "firebrick4",show_arrow = F) +
#' scale_fill_gradient2(low="dodgerblue3", mid="lightgray",high="salmon3")
#' 
plot_field_direction2 <- function(
  arrow_df,
  seurat_obj,
  mode=c("point","grid"),
  point="E_strength",
  show_arrow=TRUE,
  scale_factor=100,
  #arrow_sf=NULL,
  #line_sf=NULL,
  grid_density=NULL,
  image=c("blank","HE","cluster"),
  arrow_color=NULL,
  point_color="navy",
  arrow_normalize=FALSE,
  ...
  #alpha=c(0.5,1),
  ){
  require(ggquiver)
  #check arg
  mode <- match.arg(mode)
  image <- match.arg(image)
  #create the data.frame for drawing

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
  }else{
    #arrow_sf <- 0
    #line_sf <- Inf
    scale_factor=0
  }

  if(is.null(arrow_color)){
    if(image=="blank"){
      arrow_color <- "firebrick3"
    }else{
      arrow_color <- "black"
    }
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
  switch(image,
    blank={
      p_arrow <- 
        ggplot(data = arrow_df,aes(x=x,y=y)) +
          geom_point(
            data=arrow_df,
            aes(x=x,y=y,fill=.data[[point]]),
            size=3, shape = 23, stroke = 0.1
          )+ 
          scale_fill_gradient(low="lightgrey", high=point_color) +
          geom_quiver(
            data=arrow_draw,
            aes(x=x,y=y,u=scale_factor*Ex, v=scale_factor*Ey),
            linewidth = 1, color = arrow_color
            #arrow = arrow(length = unit(arrow_size, "npc"))
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
          geom_quiver(
            data=arrow_draw,
            aes(x=x,y=y,u=scale_factor*Ex, v=scale_factor*Ey,fill=NULL),
            linewidth = 1, color = arrow_color
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
          geom_quiver(
            data=arrow_draw,
            aes(x=x,y=y,u=scale_factor*Ex, v=scale_factor*Ey,fill=NULL),
            linewidth = 1, color = arrow_color
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
#'
#' @return return a ggplot2 object plot
#' @examples
#' LR_selected <- 
#'   c("Fgf9.Fgfr1","Fgf9.Fgfr2","Fgf9.Fgfr3",
#'     "Fgf7.Fgfr1","Fgf7.Fgfr2",
#'     "Fgf5.Fgfr2","Fgf5.Fgfr3"
#'   )
#' plot_db_score_dot(
#'   db_C2C_score_list,complex_data$kept_db,
#'   #LR_family = c("FGF"),
#'   LR_pair = LR_selected,
#'   source_use = c("7","11"),
#'   target_use = c("7","9","11"),
#'   scale = TRUE
#'   )
plot_db_score_dot <- function(
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
    filter(p_value<p_val) %>%
      mutate(
        LR_pair=str_replace(.data[["id"]],pattern = "\\.",replacement = "->"),
        S2T=paste0(.data[["Source"]],"->",.data[["Target"]])
        )
  if(scale && nrow(draw_df)!=1){
    draw_df <- draw_df %>% 
      mutate(scale_score = scale(.data[["raw_score"]]))
    score_title <- "scaled score"
  }else{
    if(nrow(draw_df)==1){message("Only one row in filtered data. Will not scale.")}
    draw_df <- draw_df %>% 
      mutate(scale_score = .data[["raw_score"]])
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
MD2_color_picker <- function(n){
  require(ggsci)
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
      preAllocateTracks = list(track.height = max(strwidth(chord_order))),
      link.arr.type = "big.arrow",
      scale = FALSE
      )

  circos.track(track.index = 2, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      xplot = get.cell.meta.data("xplot")
      ylim = get.cell.meta.data("ylim")
      sector.name = str_split(get.cell.meta.data("sector.index"),pattern = "@",simplify = TRUE)[,2]
      circos.text(mean(xlim), ylim[2], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5),cex = 1)
    }, bg.border = NA) # here set bg.border to NA is important
  title(main_title)

  highlight.sector(chord_order[str_detect(chord_order,pattern = "^S")], track.index = 1,col = NA,
      text = "Sender",cex = 1.2, facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  highlight.sector(chord_order[str_detect(chord_order,pattern = "^R")], track.index = 1,col = NA,
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
#' @examples
#' plot_LR_cluster_chord(
#'  db_C2C_score_list,complex_data$kept_db,
#'    LR_pair = c("Fgf1.Fgfr1","Fgf1.Fgfr2"),
#'    #LR_family=c("FGF","TGFb")
#'    source_use=c("1","10","11","2","9"),
#'    target_use = c("7","9","11","12")
#'  )
plot_LR_cluster_chord <- function(
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
  require(circlize)
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
      filter(p_value<0.05) %>%
      #filter(Source!=Target) %>%
      mutate(Source = paste0("S@",Source),Target = paste0("R@",Target)) 

  #print(pic_df)

  temp <- which(chord_order %in% c(pic_df$Source,pic_df$Target))
  chord_order <- chord_order[temp]
  grid_col <- grid_col[temp]
  chord_group <- str_split(chord_order,pattern = "@",simplify = TRUE)[,1]
  names(chord_group) <- chord_order

  title0 <- NULL
  if(!is.null(LR_family)){title0 <- paste0(LR_family," family signaling")}
  if(length(LR_pair)==1){
    title0 <- str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }
  if(!is.null(title)){
      title <- title
  }else{title <- title0}

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
chord_gene_core_function <- function(
  pic_df,grid_col,legend_col,chord_order,chord_group,main_title
){
  circos.clear()
  circos.par(start.degree = 176)
  chordDiagram(pic_df[,c("source_lig","target_rec","raw_score")],
      grid.col = grid_col,
      #col=link_color,
      directional = 1,big.gap = 10,
      direction.type = c("diffHeight", "arrows"),
      group = chord_group,
      order=chord_order,
      annotationTrack = "grid", 
      preAllocateTracks = list(track.height = max(strwidth(chord_order))),
      link.arr.type = "big.arrow",
      scale = FALSE
      )

  circos.track(track.index = 2, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      xplot = get.cell.meta.data("xplot")
      ylim = get.cell.meta.data("ylim")
      temp = str_split_fixed(get.cell.meta.data("sector.index"),pattern = "[@\\.]",n=3)
      sector.name = temp[,3]
      circos.text(mean(xlim), ylim[2], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(-0.1, 0.5),cex = 1)
    }, bg.border = NA) # here set bg.border to NA is important

  highlight.sector(chord_order[str_detect(chord_order,pattern = "^S")], track.index = 1,col = NA,
      text = "Sender",cex = 1.2, facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  highlight.sector(chord_order[str_detect(chord_order,pattern = "^R")], track.index = 1,col = NA,
      text = "Receiver",cex = 1.2,facing = "bending.outside", niceFacing = TRUE, text.vjust = 2.2)
  
  legend <- ComplexHeatmap::Legend(at = names(legend_col), type = "grid", legend_gp = grid::gpar(fill = legend_col), title = "Cluster")
  ComplexHeatmap::draw(legend, x = unit(1, "npc")-unit(20, "mm"), just = c("right"))

  title(main_title)
  
  circos.clear()
}



#' Chordigram plot for genes in Ligand-receptor pairs
#'
#' Chordigram plot for genes in Ligand-receptor pairs
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
#' @examples
#' plot_LR_gene_chord(
#'  db_C2C_score_list,complex_data$kept_db,
#'    LR_pair = c("Fgf1.Fgfr1","Fgf1.Fgfr2"),
#'    #LR_family=c("FGF","TGFb")
#'    source_use=c("1","10","11","2","9"),
#'    target_use = c("7","9","11","12")
#'  )
plot_LR_gene_chord <- function(
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
  require(circlize)
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
      filter(p_value<0.05) %>%
      #filter(Source!=Target) %>%
      mutate(Source = paste0("S@",Source),Target = paste0("R@",Target)) %>%
      mutate(Source=str_replace_all(Source,pattern = "\\.",replacement = "_"),
            Target=str_replace_all(Target,pattern = "\\.",replacement = "_")) %>%
      #mutate(Source = paste0(Source,"_",id),Target = paste0(Target,"_",id))
      mutate("lig"=str_split(id,pattern = "\\.",simplify = TRUE)[,1],
            "rec"=str_split(id,pattern = "\\.",simplify = TRUE)[,2]) %>%
      mutate(source_lig=paste0(Source,".",lig),
            target_rec=paste0(Target,".",rec))


  #** get the color of each cluster
  if(is.null(cluster_color)){
    clu_color <- assign_clu_col_lite(C2C_score_df)
  }else{
    clu_color <- check_cluster_color(cluster_color,C2C_score_df,source_use,target_use)
  }

  #** assign the color of each cluster and ligand/receptor
  color_df <- clu_color %>% as.data.frame()
  colnames(color_df) <- "color"
  color_df$clu <- rownames(color_df)
  LR_clu <- c(pic_df$source_lig,pic_df$target_rec) %>% table() %>% names() %>%
    str_split_fixed(pattern="[@\\.]",n=3) %>% as.data.frame()
  colnames(LR_clu) <- c("SR","clu","LR")

  LR_clu %<>% left_join(color_df,by="clu")
  grid_col <- LR_clu$color
  chord_order <- 
    names(grid_col) <- c(pic_df$source_lig,pic_df$target_rec) %>% table() %>% names()


  #** filter out not used color,and this arg is used in legend
  clu_color_use <- clu_color[LR_clu$clu %>% table %>% names]

  #** filter out not used color,and those args is used in chord plot
  temp <- which(chord_order %in% c(pic_df$source_lig,pic_df$target_rec))
  chord_order <- chord_order[temp]
  grid_col <- grid_col[temp]
  chord_group <- str_split(chord_order,pattern = "@",simplify = TRUE)[,1]
  names(chord_group) <- chord_order

  #** set title
  title0 <- NULL
  if(!is.null(LR_family)){title0 <- paste0(LR_family," family signaling")}
  if(length(LR_pair)==1){
    title0 <- str_replace(LR_pair,pattern = "\\.",replacement = "->")
  }
  if(!is.null(title)){
      title <- title
  }else{title <- title0}

  chord_gene_core_function(
      pic_df,
      grid_col,clu_color_use,
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
#' @examples
#' plot_LR_cluster_heatmap(
#'   db_C2C_score_list,
#'   complex_data$kept_db,
#'   #C2C_score_df,
#'   LR_pair = c("Fgf1.Fgfr1","Fgf1.Fgfr2"),
#'   source_use=c("1","10","11","2","9"),
#'   target_use = c("7","9","11","12"),
#'   title=NULL,
#'   method_use="strength",
#'   normalize=FALSE,
#'   ht_col="Reds"
#' )
plot_LR_cluster_heatmap <- function(
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
      title = paste0(ifelse(normalize,"Scale ",""),
      ifelse(method_use=="count","count","strength"))
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
    suffix_title <- str_replace(LR_pair,pattern = "\\.",replacement = "->")
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

