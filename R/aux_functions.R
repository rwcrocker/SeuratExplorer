# aux_functions.R
# auxillary functions needed for the package

#' Launch SeuratExplorer
#' @import shiny
#' @return In-browser Shiny Application launch
#' @export
launchSeuratExplorer <- function(){
  app = shinyApp(ui, server)
  runApp(app, launch.browser = TRUE)
}

#' Get legit annotations from metadata column names
#' @import purrr
#' @param df Metadata data.frame from Seurat object
#' @param max_unique Cut-off for total unique elements in annotations
legit_annotations = function(df, max_unique = 100){
  legit_type = map_lgl(df, .f = ~ is.character(.x) | is.factor(.x))
  legit_num = map_lgl(df, ~length(unique(.x)) < 100)
  la = names(df)[legit_type & legit_num]
  return(la)
}

#' Parse input of comma seperated genes
#' @import stringr dplyr magrittr
#' @param str string to parse
#' @param obj Seurat object to crossreference parsed genes
read_genes <- function(str, obj){
  genes = strsplit(str, split=",") %>%
    unlist() %>%
    stringr::str_trim() %>%
    as.vector() %>%
    unique()
  genes_in_obj = genes[genes %in% rownames(obj)]
  return(genes_in_obj)
}

#' Replace 'None' with NULL for meaningful selections
#' @param input string
replace_nones = function(input){
  out = input
  is_none = input == "None"
  if(is_none) out = NULL
  return(out)
}

#' VlnPlot wrapper for stacked/non-stacked plots
#' @import Seurat
#' @param obj Seurat object
#' @param features Vector of gene names
#' @param split.by Metadata column name to split violin
#' @param idents Vector of identies to plot
flexible_vln = function(obj, features, split.by=NULL, idents=NULL){
  is_multi = length(features) > 1
  plot = VlnPlot(obj, features = features,
          split.by = split.by, idents = idents,
          stack=is_multi, flip = is_multi)
  return(plot)
}


# # TODO figure out how to implement
# # Split auto-faceted Dimplot into list of individual ggplots
# split_facet = function(x){
#   facet_vars = names(x$facet$params$facets)        
#   x$facet = ggplot2::ggplot()$facet             
#   datasets = split(x$data, x$data[facet_vars])   
#   new_plots = lapply(datasets, function(new_data){
#     x$data = new_data
#     return(x)
#     })
# } 
