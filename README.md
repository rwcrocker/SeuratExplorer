# SeuratExplorer
An interactive R shiny application for exploring scRNAseq data processed in Seurat

## Installation
To use SeuratExplorer, please install the required prerequisites.
```{r}
# Install prerequisite packages
install.packages(c("Seurat", "ggplot2", "shiny", "shinydashboard", "tidyverse", "scales", "cowplot"))
```
Install SeuratExplorer from GitHub.
```{r}
# Install SeuratExplorer
install.packages('remotes')
remotes::install_github(repo = "Hla-Lab/SeuratExplorer")
```
## Running the App
To open SeuratExplorer, launch the app from the R console or Rstudio
```{r}
launchSeuratExplorer()
```
