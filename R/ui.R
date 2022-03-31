# ui.R
# R shiny UI for SeuratExplorer

#' @import Seurat ggplot2 cowplot
#' @import dplyr scales shiny shinydashboard
#' @import stringr

#Defaults
ex_genes = "Actb, Tubb1, S1pr1..."

## Header
header = dashboardHeader(title = "SeuratExplorer")

## Sidebar

sidebar = dashboardSidebar(
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


body = dashboardBody(div(class= "tab-content", tab_list))

ui = dashboardPage(header, sidebar, body)





