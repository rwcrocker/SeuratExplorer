################################################################################
#              Functions for Seurat Explorer                                   #
################################################################################


#' @import Seurat ggplot2 cowplot
#' @import dplyr scales shiny shinydashboard
#' @import stringr
NULL

### App Fxns  ###

#' Launch SeuratExplorer
#' @return In-browser Shiny Application launch
#' @export
launchSeuratExplorer <- function(){
  app = shinyApp(ui, server)
  runApp(app, launch.browser = TRUE)
}


###  Seurat Related Fxns  ###

switch_idents <- function(obj, meta.col){
  obj@active.ident <- factor(obj@meta.data[[meta.col]])
  Idents(obj) <- obj@active.ident
  return(obj)
}

read_genes <- function(str){
  genes = as.vector(stringr::str_trim(unlist(strsplit(str, split=","))))
  return(genes)
}

###  ggplot2 themes  ###

dotplot_theme <- theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text = element_text(family = "Arial", face = "bold"),
                       panel.grid.major = element_line(color = "#e0e0e0"))

bare_theme <- theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text = element_text(family = "Arial", face = "bold"))
