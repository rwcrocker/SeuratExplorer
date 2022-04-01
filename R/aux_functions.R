# aux_functions.R
# auxillary functions needed for the package

#' @import Seurat ggplot2 cowplot
#' @import dplyr scales shiny shinydashboard
#' @import stringr purrr
NULL

#' Launch SeuratExplorer
#' @return In-browser Shiny Application launch
#' @export
launchSeuratExplorer <- function(){
  app = shinyApp(ui, server)
  runApp(app, launch.browser = TRUE)
}


## Internal Fxns

#' Get legit annotations from metadata colnames
legit_annotations = function(df, max_unique = 100){
  legit_type = map_lgl(df, .f = ~ is.character(.x) | is.factor(.x))
  legit_num = map_lgl(df, ~length(unique(.x)) < 100)
  la = names(df)[legit_type & legit_num]
  return(la)
}

# Parse input of comma seperated genes
read_genes <- function(str, obj){
  genes = strsplit(str, split=",") %>%
    unlist() %>%
    stringr::str_trim() %>%
    as.vector() %>%
    unique()
  genes_in_obj = genes[genes %in% rownames(obj)]
  return(genes_in_obj)
}

# Replace None with NULL for meaningful selections
replace_nones = function(input){
  out = input
  is_none = input == "None"
  if(is_none) out = NULL
  return(out)
}

# VlnPlot wrapper wautomated for stacked/non-stacked plots
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
