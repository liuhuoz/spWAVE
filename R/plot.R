#' Plot field projection on image
#'
#' Plot given LR field projection over image of HE，cluster,
#' expression or other quantity.
#'
#' @param arrow_df data.frame,contain barcode, coordinates, field vector, expression and other quantity. 
#' Usually using the output from \code{\link{extract_LR_field_result}}
#' @param seurat_obj Seurat object, used to extract image.
#' @param mode arrow show in each point(spot) or in grid. 
#' Only support "point" or "grid", default is "point".
#' @param point point(spot) value to show, must be in colnames of arrow_df.
#' @param show_arrow logical, whether to show arrow, default is TRUE.
#' @param arrow_sf arrow head size scale factor, smaller value means larger head size.
#' @param line_sf arrow line size scale factor, smaller value means larger line size.
#' @param grid_density grid density of arrow, smaller value means sparser arrow over image,default is 0.5.
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
#' plot_field_direction(hgin_IGF,lze22_hgin,lze22_hgin,
#'  pt.size.factor=100,alpha=0.7, #pt.size.factor and alpha are passed to SpatialPlot
#'  image="cluster",mode="grid",
#'  show_arrow = T,arrow_color = "cyan2")
#'  
#' # plot field projection over cluster image with custom arrow color
#' # cluster info was extracted from \code{\link{seurat_obj@active.ident}}
#' plot_field_direction(mouse_fgfr_str,seurat_obj = brain,
#'  image="cluster",mode='grid',arrow_color="cyan2")
#' 
#' # plot field with custom point fill color indicating quantity in spot
#' plot_field_direction(temp,lze22_nor,point = "E_strength",mode="grid",
#'  line_sf = 0.3,arrow_sf = 0.05,
#'  point_color="salmon3",arrow_color = "black",show_arrow = F)
#' 
#' # plot field with custom triple colors indicating quantity in spot
#' plot_field_direction(temp,lze22_nor,point = "Rel_LR_Exp",mode="grid",
#'  line_sf = 0.3,arrow_sf = 0.05,arrow_color = "firebrick4",show_arrow = F) +
#'    scale_fill_gradient2(low="dodgerblue3", mid="lightgray",high="salmon3")
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
