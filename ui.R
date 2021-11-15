# SECexplorer MCF7x3
# MH, TDM, LR
#################################################################################
# User Interface
################

## Load packages
library(shiny)
library(ggplot2)
library(CCprofiler)
library(RColorBrewer)
library(plotly)
library(DT)

## Prepare data

# load data
for (i in list.files("data/")){
  assign(gsub(".rds", "", i), readRDS(paste0("data/",i)))
}

# assign annotation columns to vectors
annotations <- c("Protein.Group", "Genes", "Corum_complex_name", "Corum_complex_id")
Protein.Group = unique(data_prot$Protein.Group)
Genes = unique(data_prot$Genes)
Corum_complex_name = unique(corum2017$complex_name)
Corum_complex_id = unique(corum2017$complex_id)


## Define User Interface 
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("SECexplorer MCF7x3"),
  
  # Sidebar to select proteins for display
  sidebarPanel(
    selectInput("fcolumn", "Choose annotation type for protein/gene selection", annotations, selected = "Genes"),
    
    uiOutput("fcolumnvalues"), 
    #The possible choices of this field are calculated on the server side and passed over by uiOutput
    
    checkboxInput("show_monomer_markers", label = "Show monomer expected elution fraction markers", value = TRUE),
    
    textInput("legend.position", label = "Plot legend position (or \"none\")", value = "right"),
    
    p("Select genes/protein of interest above to display elution profiles."),
    p("App author: M. Heusel, moritz.heusel@med.lu.se"),
  ),
  
  # mainPanel: Show interactive line plot the selected Traces
  mainPanel(
    plotlyOutput("prot_chrom"),
    plotlyOutput("pep_chrom"),
    DT::dataTableOutput("table")
    
  )
))
