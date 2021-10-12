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
  shinyApp(ui, server)
}


###  Seurat Related Fxns  ###

switch_idents <- function(obj, meta.col){
  obj@active.ident <- factor(obj@meta.data[[meta.col]])
  Idents(obj) <- obj@active.ident
  return(obj)
}

read_genes <- function(str){
  genes = as.vector(str_trim(unlist(strsplit(str, split=","))))
  return(genes)
}


###  General Fxns  ###

check_packages_fxns = function(req_packages=NULL, req_fxns=NULL){
  packages_loaded=unlist(lapply(req_packages, require, character.only = TRUE))
  fxns_loaded=unlist(lapply(req_fxns, existsFunction))
  missing_packages=c()
  missing_fxns=c()
  if(!(all(packages_loaded) & all(fxns_loaded))){
    missing_packages = req_packages[which(packages_loaded == FALSE)]
    missing_fxns = req_fxns[which(fxns_loaded == FALSE)]
    stop(paste0("ERROR: Missing packages/fxns: ", c(missing_packages, missing_fxns)))
  }
  else
    paste0("All required packages and functions successfully loaded...")
}


###  ggplot2 themes  ###

dotplot_theme <- theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       axis.text = element_text(family = "Arial", face = "bold"),
                       panel.grid.major = element_line(color = "#e0e0e0"))

bare_theme <- theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text = element_text(family = "Arial", face = "bold"))
