# SECexplorer MCF7x3
# MH, TDM, LR
#################################################################################
# Server functions
##################

## Load packages
library(shiny)
library(ggplot2)
library(CCprofiler)
library(RColorBrewer)
library(plotly)
library(DT)

# load data
for (i in list.files("data/")){
  assign(gsub(".rds", "", i), readRDS(paste0("data/",i)))
}

# assign annotation columns to vectors
annotations <- c("Protein.Group", "Genes", "Corum_complex_name", "Corum_complex_id")
Protein.Group = unique(data_prot$Protein.Group)
Genes = unique(data_prot$Genes)
Corum_complex_name = unique(corum2017[, n_detected:=length(unique(protein_id %in% Protein.Group)), complex_name][n_detected>1]$complex_name)
Corum_complex_id = unique(corum2017[, n_detected:=length(unique(protein_id %in% Protein.Group)), complex_name][n_detected>1]$complex_id)
n_frac = max(as.numeric(data_prot$fraction_number))

# Define server fundtions
shinyServer(function(input, output) {
  
  ## Generate Reactive Filter Value Field for UI, depending on filter column chosen
  output$fcolumnvalues <- renderUI({
    values <- sort(unique(get(input$fcolumn)))
    # values <- values[nchar(values)>0]
    selectizeInput("fvalue", "Choose proteins/genes or corum complex to display (type to search)",
                   values,
                   multiple = TRUE, 
                   options = list(maxOptions = 60000),
                   selected = sort(Genes[grep("COPS", Genes)]))
  })
  
  ## generate selected protein SEC chromatogram plot
  output$prot_chrom <- renderPlotly({
    
    # Subset data to target proteins matching the selection
    if (input$fcolumn %in% "Corum_complex_name"){
      target_proteins = corum2017[complex_name %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Corum_complex_id") {
      target_proteins = corum2017[complex_id %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Genes"){
      target_proteins = data_prot[, .(Protein.Group, Genes)][Genes %in% input$fvalue, unique(Protein.Group)]
    } else {
      target_proteins = input$fvalue
    }
      
    # subset the data for these, including Mass column
    data_prot_s = data_prot[Protein.Group %in% target_proteins]
    
    # Plot
    ncolors = length(target_proteins)
    
    prot_chrom = ggplot(data_prot_s[fraction_number>0], aes(fraction_number, y = intensity, col = Genes)) +
      facet_grid(paste("replicate",replicate)~toupper(cell_line)) +
      theme_classic() +
      ggtitle("Protein/Subunit elution profiles, top3 sum") +
      scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(ncolors))   +
      theme(legend.position=input$legend.position) 
    
    if (is.null(calibration_functions)){
      prot_chrom <- prot_chrom + geom_line()
      
    } else {
      lx.frc <- seq(10,n_frac-1,10)
      lx <- paste( lx.frc , round(calibration_functions$FractionToMW(lx.frc), 1) , sep = '\n' )
      
      prot_chrom <- prot_chrom + geom_line() +
        scale_x_continuous(breaks = lx.frc , labels = lx) + xlab("fraction\n(apparent MW [kDa])")
      
      if(input$show_monomer_markers == TRUE){
        prot_chrom <- prot_chrom + geom_point(aes(x = calibration_functions$MWtoFraction(protein_mw),
          y = -10, color = Genes), pch = 23, fill = "white", size =2)
      }
      ggplotly(prot_chrom)
    }
  })
      
  ## generate selected protein SEC chromatogram plot
  output$pep_chrom <- renderPlotly({
    
    # Subset data to target proteins matching the selection
    if (input$fcolumn %in% "Corum_complex_name"){
      target_proteins = corum2017[complex_name %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Corum_complex_id") {
      target_proteins = corum2017[complex_id %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Genes"){
      target_proteins = data_prot[, .(Protein.Group, Genes)][Genes %in% input$fvalue, unique(Protein.Group)]
    } else {
      target_proteins = input$fvalue
    }
      
    # subset the data for these, including Mass column
    data_pep_s = data_pep[Protein.Group %in% target_proteins]
    
    # Plot
    ncolors = length(target_proteins)
    
    pep_chrom = ggplot(data_pep_s[fraction_number>0], aes(fraction_number, y = intensity, col = Genes, group = Modified.Sequence)) +
      facet_grid(paste("replicate",replicate)~toupper(cell_line)) +
      theme_classic() +
      ggtitle("Peptide level elution profiles, precursor intensity sum") +
      scale_color_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(ncolors))   +
      theme(legend.position=input$legend.position) 
    
    if (is.null(calibration_functions)){
      pep_chrom <- pep_chrom + geom_line()
      
    } else {
      lx.frc <- seq(10,n_frac-1,10)
      lx <- paste( lx.frc , round(calibration_functions$FractionToMW(lx.frc), 1) , sep = '\n' )
      
      pep_chrom <- pep_chrom + geom_line() +
        scale_x_continuous(breaks = lx.frc , labels = lx) + xlab("fraction\n(apparent MW [kDa])")
      
      if(input$show_monomer_markers == TRUE){
        pep_chrom <- pep_chrom + geom_point(aes(x = calibration_functions$MWtoFraction(protein_mw),
          y = -10, color = Genes), pch = 23, fill = "white", size =2)
      }
      ggplotly(pep_chrom)
    }
  })
  
  output$table <- DT::renderDataTable({
    # Collect indices of the proteins matching the selection 
    # Subset data to target proteins matching the selection
    if (input$fcolumn %in% "Corum_complex_name"){
      target_proteins = corum2017[complex_name %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Corum_complex_id") {
      target_proteins = corum2017[complex_id %in% input$fvalue, unique(protein_id)]
    } else if (input$fcolumn %in% "Genes"){
      target_proteins = data_prot[, .(Protein.Group, Genes)][Genes %in% input$fvalue, unique(Protein.Group)]
    } else {
      target_proteins = input$fvalue
    }
    protein_annotation[Entry %in% target_proteins]
  },
    options = list(
  autoWidth = TRUE))
})
