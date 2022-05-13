# server.R
## R shiny server side for SeuratExplorer

#' Server for SeuratExplorer shiny app
#' @param input Input from the UI
#' @param output Output to send back to UI
#' @import Seurat
#' @import shiny shinydashboard ggplot2
#' @import stringr purrr dplyr scales
server = function(input, output){
  options(shiny.maxRequestSize=5*1024^3)
  
  ## Dataset tab ----
  data = reactiveValues(obj = NULL, legit_annotations=NULL, annotation_idents=NULL, loaded=FALSE)
  
  observe({
    shiny::req(input$dataset_file)
    ext = tools::file_ext(input$dataset_file$datapath)
    validate(need(ext=="rds", "Please upload a .rds file"))
    data$obj = readRDS(file = input$dataset_file$datapath)
  })
  
  # Find legit_annotations
  observe({
    req(data$obj)
    data$legit_annotations = legit_annotations(data$obj@meta.data)
    data$loaded = !is.null(data$obj)
  })
  
  # Get and set annotation / idents
  observe({
    updateSelectInput(inputId = "annotation", choices = data$legit_annotations)
  })
  
  observe({
    req(input$annotation)
    data$obj = SetIdent(data$obj, value = input$annotation)
    data$annotation_idents = as.vector(unique(data$obj@active.ident))
  })
  
  # Render metadata table
  output$dataset_meta = renderDataTable({
    req(data$obj)
    data$obj@meta.data[,data$legit_annotations]
  })
  
  # Add example/reference idents
  output$annotation_idents = renderText({
    paste0(
      paste(
        head(data$annotation_idents, 3), collapse = " "
      ),
      "..."
    )
  })
  
  # Conditional panel control based on loaded obj
  output$file_loaded = reactive({
    return(data$loaded)
  })
  outputOptions(output, 'file_loaded', suspendWhenHidden=FALSE)
  
  ## UMAP tab ----
  umap = reactiveValues(genes=NULL, splitby=NULL)
  
  observe({
    updateSelectInput(inputId = "umap_splitby", choices = c("None", data$legit_annotations), selected = "None")
  })
  
  observe({
    req(input$umap_gene)
    umap$genes = read_genes(input$umap_gene, data$obj)
  })
  
  observe({
    req(input$umap_splitby)
    umap$splitby = replace_nones(input$umap_splitby)
  })
  
  output$umap = renderPlot({
    req(data$obj)
    DimPlot(data$obj, label=FALSE, split.by = umap$splitby)
  })
  
  output$featureplot = renderPlot({
    req(umap$genes)
    FeaturePlot(data$obj, features = umap$genes, split.by = umap$splitby, combine = TRUE)
  })
  
  ## VLN tab ----
  vln = reactiveValues(genes=NULL, splitby=NULL)
  
  observe({
    updateSelectInput(inputId = "vln_splitby", choices = c("None", data$legit_annotations), selected = "None")
  })
  
  observe({
    req(input$vln_splitby)
    vln$splitby = replace_nones(input$vln_splitby)
  })
  
  observe({
    req(input$vln_gene)
    vln$genes = read_genes(input$vln_gene, data$obj)
  })
  
  observe({
    updateSelectInput(inputId = "vln_idents", choices = data$annotation_idents)
  })
  
  output$vln = renderPlot({
    req(vln$genes, data$obj)
    flexible_vln(obj = data$obj,
                 features = vln$genes,
                 split.by = vln$splitby,
                 idents = input$vln_idents
    )
  })
  
  ## DotPlot tab ----
  dot = reactiveValues(genes=NULL)
  
  observe({
    updateSelectInput(inputId = "dot_idents", choices = data$annotation_idents)
  })
  
  observe({
    req(input$dot_gene, data$obj)
    dot$genes = read_genes(input$dot_gene, data$obj)
  })
  
  output$dot = renderPlot({
    req(dot$genes, data$obj)
    DotPlot(data$obj,
            features = dot$genes,
            cols = c("white", "black"),
            idents = input$dot_idents)+
      RotatedAxis()
  })
  
  ## DE tab ----
  de = reactiveValues(conditions=NULL, ref_cells=NULL, exp_cells=NULL, table=NULL, loaded = NULL)
  
  # Update choices for annotation and DE varible
  observe({
    updateSelectInput(inputId = "de_ident", choices = data$annotation_idents)
  })
  observe({
    updateSelectInput(inputId = "de_variable", choices = data$legit_annotations)
  })
  
  
  # Update choices for condition
  observe({
    req(data$obj, input$de_variable)
    de$conditions = data$obj@meta.data[input$de_variable] %>%
      unique() %>%
      as.vector()
  })
  observe({
    updateSelectInput(inputId = "de_reference",
                      choices = de$conditions)
  })
  observe({
    updateSelectInput(inputId = "de_experimental",
                      choices = de$conditions[de$conditions != input$de_reference])
  })
  
  # Warning
  output$de_warning = renderText({
    paste0('Differential Expression testing is computationally intensive and may take some time. Only press the "Analyze" button once!')
  })
  
  observeEvent(input$de_analyze, {
    de$table = FindMarkers(object = data$obj,
                           ident.1 = input$de_experimental,
                           ident.2 = input$de_reference,
                           group.by = input$de_variable,
                           subset.ident = input$de_ident,
                           only.pos = FALSE)
    de$table = cbind(rownames(de$table), de$table)
    colnames(de$table)[1] = "gene_name"
    de$loaded = ~is.null(de$table)
  })
  
  output$de_table = renderDataTable({
    de$table
  })
  
  
  # Control conditional download button
  output$de_loaded = reactive({
    return(de$loaded)
  })
  outputOptions(output, 'de_loaded', suspendWhenHidden=FALSE)
  
  # Handle downloads of DEGs as csv
  output$de_download = downloadHandler(
    filename = function() {
      paste(input$de_reference, "_vs_",input$de_experimental, "_in_", input$de_ident, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(de$table, file = file, row.names = FALSE)
    }
  )
  
  ## LR tab ----
  
  lr = reactiveValues(table=NULL, subobj=NULL, condition=NULL, conditions=NULL, signatures=NULL, degs=NULL, db=NULL)
  
  output$lr_warning = renderText({
    paste0('Ligand-Receptor analysis is computationally intensive and may take some time. Be conservative. Only press the "Analyze" button once!')
  })
  
  observe({
    updateSelectInput(inputId = "lr_db", choices = c("Mouse", "Human"), selected = "Mouse")
  })
  
  observe({
    updateSelectInput(inputId = "lr_idents", choices = data$annotation_idents)
  })
  
  observe({
    updateSelectInput(inputId = "lr_condition", choices = c("None", data$legit_annotations), selected = "None")
  })
  
  observe({
    req(input$lr_condition)
    lr$condition = replace_nones(input$lr_condition)
  })
  
  observe({
    req(lr$condition)
    lr$conditions = data$obj@meta.data %>%
      pull(input$lr_condition) %>%
      unique() %>%
      as.vector()
    updateSelectInput(inputId = "lr_reference", choices = lr$conditions)
  })
  
  observe({
    updateSelectInput(inputId = "lr_experimental",
                      choices = lr$conditions[lr$conditions != input$lr_reference])
  })
  
  observe({
    lr$db = switch(input$lr_db,
                   Mouse = SeuratExplorer:::Mouse_CellTalkDB,
                   Human = SeuratExplorer:::Human_CellTalkDB
    )
  })
  
  observeEvent(input$lr_analyze, {
    lr$subobj = subset(data$obj, idents = input$lr_idents)
    lr$signatures = FindAllMarkers(lr$subobj, only.pos = TRUE)
    lr$table = analyze_LR(lr$signatures, lr$db)
    
    if (input$lr_condition != "None"){
      lr$degs = get_all_degs(obj = lr$subobj,
                             annotation = input$annotation,
                             condition = input$lr_condition,
                             reference = input$lr_reference,
                             experimental = input$lr_experimental,
                             padj.cut=input$lr_deg_cutoff,
                             log2fc.cut=0.1)
      lr$table = crossreference_degs(lr_table = lr$table, deg_table = lr$degs, db = lr$db)
      }
    })
  
  output$lr_table = renderDataTable(
    lr$table,
    options = list(scrollX = TRUE)
  )
  
  output$lr_download = downloadHandler(
    filename = function() {
      paste0(paste(input$lr_idents, collapse = "_"), "_by_", input$lr_condition, "_LR_table", ".csv")
    },
    content = function(file) {
      write.csv(de$table, file = file, row.names = FALSE)
    }
  )
  
  output$lr_loaded = reactive({
    return(!is.null(lr$table))
  })
  outputOptions(output, 'lr_loaded', suspendWhenHidden=FALSE)
  
  
  ## DEV ----
  output$dev_out = renderText({
    "test"
  })
  
}


