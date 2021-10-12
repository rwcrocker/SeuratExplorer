################################################################################
###                        Seurat-Explorer App                               ###
################################################################################


req_packages = c("Seurat", "ggplot2", "shiny", "shinydashboard", "dplyr", "scales", "cowplot", "stringr")
lapply(req_packages, require, character.only = TRUE)


#source("./aux_functions.R")
#source("./ligand_receptor.R")

#Defaults
ex_genes = "Actb, Tubb1, S1pr1..."

################################################################################
###                                   UI                                     ###

## Header  ###

header <- dashboardHeader(title = "SeuratExplorer")

###  Sidebar  ###

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Dataset", tabName = "dataset", icon=icon(name="database"))
    ),
  conditionalPanel(
    condition = "output.file_loaded",                                         
    selectInput("annotation", label = "Annotation:", choices = NULL),
    textOutput("annotation_identities")
  ),
  sidebarMenu(
    menuItem("UMAP", tabName = "umap", icon=icon(name="cloudsmith")),
    menuItem("Violin", tabName = "vln", icon=icon("hourglass")),
    menuItem("DotPlot", tabName = 'dot', icon=icon("braille")),
    menuItem("Differential Expression", tabName = 'de', icon=icon("list")),
    menuItem("Ligand-Receptor", tabName = 'lr', icon=icon("connectdevelop")),
    menuItem("Developer Testing", tabName = 'dev', icon=icon("braille"))
  )
)

###  BODY  ###

# Pages

tab_list=list()
tab_list[["dataset"]]=
  tabItem("dataset",
          box(
            fileInput("dataset_file", "Choose Seurat .rds file:", accept = '.rds'), width = 4
            ),
          dataTableOutput("dataset_meta")
          )

tab_list[["umap"]]=
  tabItem("umap",
          sidebarLayout(
            position = "left", fluid = TRUE,
            mainPanel(
                plotOutput("umap"),
               plotOutput("featureplot"),
               width = 8
            ),
          sidebarPanel(
            selectInput("umap_splitby", label = "Split by:", choices = NULL),
            textInput("umap_gene", "Feature(s):", placeholder = ex_genes, value = NULL),
            width = 4
            )
          )
        )

tab_list[["vln"]]=
  tabItem("vln",
          sidebarLayout(
            position = "left", fluid = TRUE,
            mainPanel(
              plotOutput("vln")
              ),
            sidebarPanel(
              textInput("vln_gene", label = "Feature(s):", placeholder = ex_genes, value = NULL),
              selectInput("vln_splitby", label = "Split by:", choices = NULL),
              selectInput("vln_idents", label = "Identities:", choices = NULL, multiple = TRUE)
              )
            )
          )

tab_list[["dot"]]=
  tabItem("dot",
          sidebarLayout(
            position = "left", fluid = TRUE,
            mainPanel(
              plotOutput("dot")
            ),
            sidebarPanel(
              textInput("dot_gene", "Feature(s):", placeholder = ex_genes, value = NULL),
              selectInput("dot_idents", label = "Identities:", choices = NULL, multiple = TRUE)
            )
          )
        )

tab_list[["de"]]=
  tabItem("de",
          box(
            textOutput("de_warning"), title = "WARNING", background = "red", width = 12
            ),
          box(
            selectInput("de_clusters", label = "Cluster(s):", choices = NULL, multiple = FALSE),
            selectInput("de_comparison", label = "Conditions:", choices = NULL),
            selectInput("de_reference", label = "Reference condition:", choices = NULL),
            selectInput("de_experimental", label = "Experimental condition:", choices = NULL)
          ),
          box(
            actionButton("de_analyze", "Analyze"),
            width = 6),
          dataTableOutput("de_table"),
          conditionalPanel(
            condition = "output.de_loaded",                                 #Change to a confirmation of loaded de table
            downloadButton("de_download", "Download"),
            width = 4
            )
          )

tab_list[["lr"]]=
  tabItem("lr",
          box(
            textOutput("lr_warning"), title = "WARNING", background = "red", width = 12
          ),
          box(
            selectInput("lr_db", label = "Database:", choices = NULL, multiple = FALSE),
            selectInput("lr_clusters", label = "Cluster(s):", choices = NULL, multiple = TRUE),
            selectInput("lr_condition", label = "Condition (Optional):", choices = NULL, multiple = FALSE),
            conditionalPanel(
              condition = "input.lr_condition != 'None'",
              selectInput("lr_reference", label = "Reference condition:", choices = NULL),
              selectInput("lr_experimental", label = "Experimental condition:", choices = NULL)
            )
          ),
          box(  
            sliderInput("lr_signature_cutoff", label = "Signature Padj Cutoff:", min = 0.01, max = 0.15, value = 0.05, step = 0.01),
            sliderInput("lr_deg_cutoff", label = "DEG Padj Cutoff:", min = 0.01, max = 0.15, value = 0.05, step = 0.01),
          ),
          box(
            actionButton("lr_analyze", "Analyze"),
            width = 4),
          dataTableOutput("lr_table"),
          conditionalPanel(
            condition = "output.lr_loaded",
            downloadButton("lr_download", "Download"),
            width = 4
          )
        )

tab_list[["dev"]]=
  tabItem("dev",
          textInput(inputId = "test", label = "Input test")
          )


body <- dashboardBody(div(class= "tab-content", tab_list))

ui <- dashboardPage(header, sidebar, body)         


################################################################################
##                               Server                                       ##

# Defaults
categorical_col_exceptions = c("cells")                                         #exceptions to prevent crash from too many identities
lr_db_path = c("./resources/")





server = function(input, output){
    options(shiny.maxRequestSize=5*1024^3)                                      #required for large seurat objs
  
  #Initalize reactive values to hold changing data for each tab
  data <- reactiveValues(obj = NULL, meta=NULL, categorical_cols=NULL, annotation_identities=NULL, is_loaded=FALSE)
  umap <- reactiveValues(genelist=NULL, splitby=NULL, featplot=NULL)
  vln <- reactiveValues(genelist=NULL, identlist=NULL, is_multi=NULL, splitby=NULL)
  dot <- reactiveValues(genelist=NULL, identlist=NULL)
  de <- reactiveValues(comparison=NULL, ref_cells=NULL, exp_cells=NULL, table=NULL, is_loaded = FALSE)
  lr <- reactiveValues(table=NULL, subobj=NULL, condition_options=NULL, signatures=NULL, degs=NULL, is_loaded = FALSE)
  
  ###   Dataset tab   ###
  
  #Always fetch the new obj when loaded
  observe({
    shiny::req(input$dataset_file)
    ext <- tools::file_ext(input$dataset_file$datapath)
    validate(need(ext=="rds", "Please upload a .rds file"))
    data$obj = readRDS(file = input$dataset_file$datapath)
    data$is_loaded=TRUE
  })
  
  observe({
    req(data$obj)
    data$meta = data$obj@meta.data
    data$categorical_cols = colnames(data$meta)[which(sapply(data$meta, class) %in% c("character", "factor") & !(colnames(data$meta) %in% categorical_col_exceptions))]
  })
  
  
  #set categorical_cols and options for active ident
  observe({
    updateSelectInput(inputId = "annotation", choices = data$categorical_cols)
  })

  #set obj active ident
  observe({
    req(input$annotation)
    data$obj = switch_idents(data$obj, input$annotation)
    data$annotation_identities = as.vector(unique(data$obj@active.ident))
  })
  
  #render meta table
  output$dataset_meta <- renderDataTable({
    data$meta[,data$categorical_cols]
  })
  
  #add example idents
  output$annotation_identities <- renderText({
    paste0(
      paste(
        head(data$annotation_identities, 3),collapse = " "
      ),
      "..."
    )
  })
  
  #check if a dataset is loaded, control conditional tab panel
  output$file_loaded = reactive({
    return(data$is_loaded)
  })
  outputOptions(output, 'file_loaded', suspendWhenHidden=FALSE)

  
###   UMAP tab    ###
  
  observe({
    updateSelectInput(inputId = "umap_splitby", choices = c("None", data$categorical_cols))
  })

  observe({
    umap$genelist = read_genes(input$umap_gene)
    umap$splitby = input$umap_splitby
    if (input$umap_splitby == "None"){
      umap$splitby = NULL
    }
  })
  
  output$umap <- renderPlot({
    req(data$obj)
    DimPlot(data$obj, label=FALSE, split.by = umap$splitby)
  })

  output$featureplot <- renderPlot({
    req(umap$genelist)
    FeaturePlot(data$obj, features = umap$genelist, split.by = umap$splitby, combine = TRUE)
  })


###   VLN tab   ###
  
  observe({
    updateSelectInput(inputId = "vln_splitby", choices = c("None", data$categorical_cols))
  })
  
  observe({
    vln$genelist = as.vector(stringr::str_trim(unlist(strsplit(input$vln_gene, split=","))))
    vln$is_multi = length(vln$genelist) > 1
    updateSelectInput(inputId = "vln_idents", choices = data$annotation_identities)
    
    vln$splitby=input$vln_splitby
    if (input$vln_splitby == "None"){
      vln$splitby = NULL
    }
  })
  
  
  output$vln <- renderPlot({
    req(vln$genelist)
    VlnPlot(data$obj, features = vln$genelist,
            stack = vln$is_multi, flip = vln$is_multi,
            split.by = vln$splitby,
            idents = input$vln_idents)+
      bare_theme
  })
  
###   DotPlot tab   ###
  
  observe({
    dot$genelist = as.vector(stringr::str_trim(unlist(strsplit(input$dot_gene, split=","))))
    updateSelectInput(inputId = "dot_idents", choices = data$annotation_identities)
  })
  
  output$dot <- renderPlot({
    req(dot$genelist)
    DotPlot(data$obj, features = dot$genelist,
                     cols = c("white", "black"), idents = input$dot_idents)+RotatedAxis()+dotplot_theme
  })
  
###   DE tab    ###
  
  observe({
    updateSelectInput(inputId = "de_clusters", choices = data$annotation_identities)
  })
  
  observe({
    updateSelectInput(inputId = "de_comparison", choices = data$categorical_cols)
  })
  
  observe({
    de$comparison = as.vector(unique(data$meta[[input$de_comparison]])) 
    updateSelectInput(inputId = "de_reference", choices = de$comparison)
    updateSelectInput(inputId = "de_experimental", choices = de$comparison)
  })
  
  output$de_warning <- renderText({
    paste0('Differential Expression testing is computationally intensive and may take some time. Only press the "Analyze" button once!')
  })
  
  observeEvent(input$de_analyze, {
    de$table = FindMarkers(object = data$obj, ident.1 = input$de_experimental, ident.2 = input$de_reference, group.by = input$de_comparison, subset.ident = input$de_clusters, only.pos = FALSE)
    de$table = cbind(rownames(de$table), de$table)
    colnames(de$table)[1] = "gene_name"
    de$is_loaded=TRUE
  })
  
  output$de_table <- renderDataTable({
    de$table
  })
  
  output$de_download <- downloadHandler(
    filename = function() {
      paste(input$de_reference, "_vs_",input$de_experimental, "_in_", input$de_clusters, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(de$table, file = file, row.names = FALSE)
    }
  )
  
  #check if a table is generated, control conditional button
  output$de_loaded = reactive({
    return(de$is_loaded)
  })
  outputOptions(output, 'de_loaded', suspendWhenHidden=FALSE)

###  LR tab  ###
  
  output$lr_warning <- renderText({
    paste0('Ligand-Receptor analysis is computationally intensive and may take some time. Be conservative. Only press the "Analyze" button once!')
  })
  
  observe({
    updateSelectInput(inputId = "lr_db", choices = dir("./resources/CellTalkDB/"))
  })
  
  observe({
    updateSelectInput(inputId = "lr_clusters", choices = data$annotation_identities)
  })
  
  observe({
    updateSelectInput(inputId = "lr_condition", choices = c("None", data$categorical_cols), selected = "None")
  })
  
  observe({
    lr$condition_options = as.vector(unique(data$meta[[input$lr_condition]])) 
    updateSelectInput(inputId = "lr_reference", choices = lr$condition_options)
    updateSelectInput(inputId = "lr_experimental", choices = lr$condition_options)
  })
  
  observeEvent(input$lr_analyze, {
    lr$subobj = subset(data$obj, idents = input$lr_clusters)
    lr$signatures = FindAllMarkers(lr$subobj, only.pos = TRUE)%>%
      filter(p_val_adj < input$lr_signature_cutoff)
    lr$table = Analyze_LR(lr$signatures, db_path = paste0("./resources/CellTalkDB/", input$lr_db))
    lr$is_loaded = TRUE
    if (input$lr_condition != "None"){
      lr$degs = Find_Condition_DEGs(lr$subobj, condition = input$lr_condition, reference.cond = input$lr_reference, experimental.cond = input$lr_experimental, padj.cut = input$lr_deg_cutoff)
      lr$table = Crossreference_LR(lr$table, DEG.table = lr$degs, db.path = paste0("./resources/CellTalkDB/", input$lr_db))
    }
  })
  
  output$lr_table <- renderDataTable(
    lr$table,
    options = list(scrollX = TRUE)
  )
  
  output$lr_download <- downloadHandler(
    filename = function() {
      paste0(paste(input$lr_clusters, collapse = "_"), "_by_", input$lr_condition, "_LR_table", ".csv")
    },
    content = function(file) {
      write.csv(de$table, file = file, row.names = FALSE)
    }
  )
  
  output$lr_loaded = reactive({
    return(lr$is_loaded)
  })
  outputOptions(output, 'lr_loaded', suspendWhenHidden=FALSE)
  
}












