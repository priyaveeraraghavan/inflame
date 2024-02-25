#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# TODOS:
# for single gene input, use autofill to help with gene names
# add a multi-gene-input page
# --> input csv
# --> volcanos --> with FDR lines (0.05), FC log2FC = +/-0.5
# --> violins comparing
# --> heatmap
# --> tab on heatmap to give different normalization (z-score or range normalized?)
# --> make sure heatmaps are separate by IP experiment 
# apex
# group --> output (hotelling) and output violin over time [try fold change vs time 0]
# single gene --> hotelling, abundance over timecourse (replicates as dots, line graph connecting [try fold change vs time 0])


# functions to write
# plot_boxplot 


library(DT)
library(ggrepel)
library(shiny)
library(tidyverse)
library(shinyjs)
library(shinydashboard)
library(pheatmap)
# Setup data

project_dir <- "/Users/priyaveeraraghavan/Documents/GitHub/inflame/Proteomics_Visualization"
load(file.path(project_dir, "data/processed_proteomics_data.RData"))
# metadata_file <- file.path(project_dir, "Proteomics_Visualization/data/designNoCTR_msSTATS.csv")
# input_data_file <- file.path(project_dir, "Proteomics_Visualization/data/all_stats_results_Lyso_1peptide.csv")

accent_colors <- c("red", "blue")
# 
# load_metadata <- function(meta_file) {
#   meta_df <- read.csv(meta_file) %>%
#     subset(condition != "Empty" | !is.na(experiment_type)) %>%
#     rowwise() %>%
#     mutate(sample_name=sprintf("%s_%s_%s", experiment_type, condition, replicate))
#   
#   return(meta_df)
# }
# meta_df <- load_metadata(metadata_file)
# 
# load_data <- function(input_file) {
#   
#   input_df <- read.csv(input_file)
#   
#   # filter out the "#" 
#   input_df <- input_df %>% 
#     subset(!grepl("#", reference) & nchar(gene_symbol) > 0) 
#   
#   # add negative log10 q values
#   for (qval in colnames(input_df)[grepl("q.val", colnames(input_df))]) {
#     input_df[,sprintf("negLog10%s", qval)] <- -log10(input_df[,qval])
#   }
#   
#   # just for now with only lyso data
#   gene_symbols <- input_df$gene_symbol
#   references <- input_df$reference
#   input_df$gene_symbol <- NULL
#   input_df$reference <- NULL
#   colnames(input_df) <- sprintf("lysoIP_%s", colnames(input_df))
#   input_df$gene_symbol <- gene_symbols
#   input_df$reference <- references
#   return(input_df)
#   
# }
# 
# input_df <- load_data(input_data_file)

# get a list of comparisons based on the metadata and on what measurements
# exist in the input data frame
# searches for log2FC columns
# get_comparisons <- function(meta_df, input_df) {
#   poss_comparisons <- data.frame(
#     expand.grid(list(condition1=unique(meta_df$condition), 
#                      condition2=unique(meta_df$condition),
#                      experiment_type=unique(meta_df$experiment_type))),
#     stringsAsFactors = FALSE) %>%
#     rowwise() %>% 
#     mutate(comparison_name=sprintf("%s_%s", 
#                                    as.character(condition1),
#                                    as.character(condition2))) %>%
#     rowwise() %>%
#     mutate(log2FC_comparison_name=sprintf("%s_log2FC_%s", 
#                                           as.character(experiment_type),
#                                           comparison_name))
#   
#   # subset for measured columns
#   valid_comparisons <- poss_comparisons %>%
#     subset(log2FC_comparison_name %in% colnames(input_df))
#   
#   return(valid_comparisons)
#   
# }
# 
# comparison_df <- get_comparisons(meta_df, input_df)


## UI-side plot definitions
define_static_plot <- function(panel, plot_type, compartment, treatvsctrl) {
  
  return(box(title=sprintf("%s %s", compartment, treatvsctrl), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_%s_static_%s_%s", panel, plot_type, compartment, treatvsctrl), 
                    width = 250, height=250)))
}

define_static_scatter <- function(panel, compartment, treatvsctrl) {
  return(box(title=sprintf("%s %s", compartment, treatvsctrl), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_scatterplot_%s_%s", panel, compartment, treatvsctrl), 
                    width = 250, height=250)))
}

define_interactive_scatter <- function(panel, compartment, treatvsctrl) {
  return(box(title=sprintf("%s %s", compartment, treatvsctrl), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_scatterplot_%s_%s", panel, compartment, treatvsctrl), 
                    width = 250, height=250, 
                    brush = brushOpts(id=sprintf("%s_plot_brush", panel), resetOnNew = TRUE))))
}


# ## Helper data wrangling 
# get_intensities_gene_group <- function(compartment_name, gene_group_df) {
#   design_subs <- subset(processed_data[['design']], compartment == compartment_name)
#   
#   # design_cols remove mean
#   design_cols <- design_subs$column_name[!grepl("Mean", design_cols$column_name)]
#   
#   gene_subset_df <- data.frame(processed_data[[compartment_name]] %>%
#                                  inner_join(gene_group_df, on='reference'))
#   
# }

get_wide_data <- function(gene_group_df, msmt_cols_to_keep) {
  # function that gets the relevant information for organelles and WCP for a gene set
  # displays in wide format
  
  cols_to_keep <- c("reference", "gene_symbol", "Annotation", "compartment", msmt_cols_to_keep)
  
  compartments <- unique(processed_data[['design']]$compartment)
  
  long_df <- do.call(rbind, 
                     lapply(compartments,
                            function(comp) inner_join(processed_data[[comp]], gene_group_df)[,cols_to_keep]))
  
  wide_df <- long_df %>% pivot_wider(id_cols=c(reference, gene_symbol, Annotation),
                                     names_from=compartment,
                                     values_from=msmt_cols_to_keep)
  
  return(wide_df)

}          
plot_heatmap <- function(compartment_name, gene_group_df) {
  # compartment_name: string compartment 
  # gene_group_df: df of ['reference', 'gene_group']
  print("Plotting heatmap")
  
  design_subs <- subset(processed_data[['design']], 
                        !grepl("Mean", column_name) & compartment == compartment_name)
  
  # design_cols remove "intensity"
  design_col_rename <- gsub("_Intensity", "", design_subs$column_name)
  
  gene_subset_df <- data.frame(processed_data[[compartment_name]] %>%
                                 inner_join(gene_group_df))
  
  heatmap_df <- gene_subset_df[,  design_subs$column_name]
  
  
  # normalize as prop of max
  #rowmaxs <- apply(heatmap_df, 1, max)
  heatmap_df <- data.frame(sweep(heatmap_df, MARGIN = 1, FUN = "/", STATS = rowMeans(heatmap_df)))
  
  # make sure rows and cols named
  rownames(heatmap_df) <- gene_subset_df$reference
  colnames(heatmap_df) <- design_col_rename
  
  # define gene group for annotation
  annot_row <- data.frame(gene_group=gene_subset_df$gene_group)
  rownames(annot_row) <- gene_subset_df$reference
  
  # define sample group for annotation
  annot_col <- data.frame(design_subs[,c("condition", "compartment")])
  rownames(annot_col) <- design_col_rename
  
  pheatmap(heatmap_df,
           annotation_row = annot_row,
           annotation_col = annot_col)
  
}

plot_scatter <- function(compartment_name, comparison, 
                         gene_group_df, annotate_genes=TRUE) {
  xcol <- paste0(c("log2FC", comparison), collapse="_")
  ycol <- paste0(c("negLog10qval", comparison), collapse="_")

  
  data_full <- data.frame(processed_data[[compartment_name]])
  data_subset <- data.frame(processed_data[[compartment_name]] %>%
    inner_join(gene_group_df))
  
  
  plt <- ggplot(data_full, aes_string(xcol, ycol)) +
    geom_point(show.legend=FALSE)  +
    geom_point(data=data_subset, 
               mapping=aes_string(xcol, ycol, color="gene_group"),
               size=2,
               show.legend = FALSE) +
    theme_classic() +
    geom_vline(xintercept=0.5, linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed")
  
  if (annotate_genes) {
    plt <- plt + geom_label_repel(data=data_subset, 
                     mapping=aes_string(xcol, ycol, color="gene_group",
                                        label="gene_symbol"),
                     size=3,
                     show.legend = FALSE)
  }
  return(plt)
}

# for gene symbols
plot_boxplot <- function(compartment_name, gene_group_df) {
  
  design_subset <- subset(processed_data[['design']], 
                          compartment == compartment_name & replicate != "Mean")
  
  column_names <- design_subset$column_name
  
  data_subset <- data.frame(processed_data[[compartment_name]] %>% 
                              inner_join(gene_group_df))
  
  if (nrow(data_subset) > 0) {

    annot_df_subs <- data_subset[,c("gene_symbol", "reference", column_names)] %>%
      pivot_longer(contains("Intensity"), 
                 names_to = "column_name", values_to="intensity") %>%
      inner_join(design_subset, by="column_name")
  
    # get the identifier for the isoform
    annot_df_subs <- annot_df_subs %>%
      rowwise() %>% mutate(identifier=strsplit(reference, "\\|")[[1]][2])
  
    # plot
    ggplot(annot_df_subs, aes(condition, intensity, color=identifier)) + 
      geom_boxplot() + 
      geom_point(position = position_dodge(width=0.8)) +
      theme_classic() +
      theme(legend.position="bottom")
  } else {
    ggplot()
  }
  
}

# for gene group 
# violin plot log2FC
plot_violin_log2FC <- function(compartment_name, gene_group_df) { #experiment, annot_df) {

  
  data_subset <- data.frame(processed_data[[compartment_name]] %>%
                              inner_join(gene_group_df))
   
  annot_df_subs <- data_subset %>% 
    dplyr::select(contains("log2FC") | matches("gene_group")) %>%
    pivot_longer(contains("log2FC"), names_to="comparison", values_to="log2FC")
  
  # modify the x axis names so they fit
  annot_df_subs <- annot_df_subs %>%
    rowwise() %>%
    mutate(comparison=gsub("vs", "\nvs\n",gsub("log2FC_", "", comparison)))
  
  # plot
  ggplot(annot_df_subs,
         aes(comparison, log2FC, color=gene_group)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
    geom_hline(yintercept=0, color='black', linetype='dashed') + 
    theme_classic() +
    theme(legend.position="bottom") +
    ggtitle(compartment_name)
}

# Define UI for application 
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Single Gene Viewer",
                              fluidRow(
                                box(width=4, textInput(inputId = "symbol", label = "Gene Symbol", 
                                                           value = "", placeholder = "NLRP3")     # might also want to try textAreaInput for longer
                                )),
                              fluidRow(
                                lapply(unique(processed_data[['design']]$compartment),
                                       function(compart) column(width=2,
                                       lapply(subset(processed_data[['comparisons']], compartment == compart)$comparison,
                                              function(comparison) { define_static_scatter("panelsymbol", compart, comparison)})
                                                               )
                                      )
                                ),
                              fluidRow(
                                lapply(unique(processed_data[['design']]$compartment),
                                       function(comp) column(width=2, define_static_plot("panelsymbol", "boxplot", comp,  "all")))
                                
                                ), 
                              
                              fluidRow(
                                box(width = 12, DT::dataTableOutput("table_symbol")))
                      ),
    tabPanel("Gene Region Explorer",
             fluidRow(
               lapply(unique(processed_data[['design']]$compartment),
                      function(compart) column(width=2,
                                               lapply(subset(processed_data[['comparisons']], compartment == compart)$comparison,
                                                      function(comparison) { define_interactive_scatter("panelregion", compart, comparison)})
                      )
               )
             ),
             #   column(width=4,
             #          lapply(subset(processed_data[['comparisons']], compartment == "lysosome")$comparison,
             #                 function(comp) { define_interactive_scatter("panelregion", "lysosome", comp)})
             #   )
             # ),
             lapply(unique(processed_data[['design']]$compartment),
                    function(comp) fluidRow(box("Gene Abundance Heatmap", 
                                                width=12, 
                                                plotOutput(sprintf("panelregion_heatmap_static_%s_all", comp))))
                    ),
             #column(width=12, define_static_plot("region", "heatmap", "lysosome", "all"))),
             #lapply(unique(processed_data[['design']]$compartment),
             #        function(comp) fluidRow(define_static_plot("region", "heatmap", comp, "all"))
             #        ),
             #  box("Gene Abundance Heatmap", width=12, plotOutput("heatmap_region"))
             #),
             fluidRow(
               box(width = 12, DT::dataTableOutput("table_region"))
             )
    ),
    tabPanel("Gene Group Viewer",
             fluidRow(
               fileInput("gene_group_input", "Gene Groupings", 
                         accept=".csv",
                         buttonLabel = "Upload")
             ),
             fluidRow(
               lapply(unique(processed_data[['design']]$compartment),
                      function(compart) column(width=2,
                                               lapply(subset(processed_data[['comparisons']], compartment == compart)$comparison,
                                                      function(comparison) { define_static_scatter("panelgroup", compart, comparison)})
                      )
               )
             ),
             fluidRow(
               lapply(unique(processed_data[['design']]$compartment),
                      function(comp) column(width=2, define_static_plot("panelgroup", "violinplot", comp,  "all")))
               
             ), 
             lapply(unique(processed_data[['design']]$compartment),
                    function(comp) fluidRow(box("Gene Abundance Heatmap", 
                                                width=12, 
                                                plotOutput(sprintf("panelgroup_heatmap_static_%s_all", comp))))
             ),
             fluidRow(
               box(width = 12, DT::dataTableOutput("table_group"))
             )
    )
  )
)

                

server <- function(input, output) {
  
  empty_df <- processed_data[['gene_annotations']] %>% 
    subset(gene_symbol == "DoesNotExist") %>%
    rowwise() %>%
    mutate(gene_group="None")

  annotation_dfs <- reactiveValues(  
    
    # gene group page
    gene_group_df = empty_df,
    # single gene page
    gene_symbol_df = empty_df,
    # gene region explorer
    gene_region_df = empty_df
  )
    
  
  observeEvent(input$symbol, ignoreNULL = TRUE, {
    if (nchar(input$symbol) > 0) {
      gene_symbols <- c(gsub(" ", "", input$symbol))
      annotation_dfs$gene_symbol_df <- subset(processed_data[['gene_annotations']], 
                                              gene_symbol %in% gene_symbols) %>%
        rowwise() %>% mutate(gene_group = "Gene_Selection") %>%
        dplyr::select(reference, gene_symbol, gene_group)
    }
    
  })
  
  observeEvent(input$panelregion_plot_brush, ignoreNULL = TRUE, {
    brush_output_id <- input$panelregion_plot_brush$outputId
    brush_compare <- gsub("panelregion_scatterplot_", "", brush_output_id)
    compartment_name <- strsplit(brush_compare, "_")[[1]][1]
    comparison_name <- strsplit(brush_compare, "_")[[1]][2]
    
    brush_xcol <- paste0(c("log2FC", comparison_name), collapse="_")
    brush_ycol <- paste0(c("negLog10qval", comparison_name), collapse="_")
    
    annotation_dfs$gene_region_df <- processed_data[[compartment_name]] %>% 
      dplyr::filter(.data[[brush_xcol]] < input$panelregion_plot_brush$xmax &
                      .data[[brush_xcol]] > input$panelregion_plot_brush$xmin &
                      .data[[brush_ycol]] < input$panelregion_plot_brush$ymax & 
                      .data[[brush_ycol]] > input$panelregion_plot_brush$ymin) %>%
      rowwise() %>% 
      mutate(gene_group="Regional_Selection") %>% 
      dplyr::select(reference, gene_symbol, gene_group)
  })
  
  observeEvent(input$gene_group_input, ignoreNULL=TRUE, {
    
    file <- input$gene_group_input
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    # read the input
    gene_gp_df <- read.csv(file$datapath)
    colnames(gene_gp_df) <- c("gene_symbol", "gene_group")
    
    # identify missing # and pct mapped, will need to display as text
    missing_ids <- setdiff(gene_gp_df$gene_symbol, processed_data[['gene_annotations']]$gene_symbol)
    num_unique_ids <- length(unique(gene_gp_df$gene_symbol))
    map_pct <- (num_unique_ids-length(missing_ids))/length(num_unique_ids)*100
    
    # merge with input_df
    annotation_dfs$gene_group_df <- processed_data[['gene_annotations']] %>% 
      inner_join(gene_gp_df)
    
  })
  
  plot_defs <- expand.grid(unique(processed_data[['design']]$compartment),
                              unique(processed_data[['comparisons']]$comparison),
                              c("symbol", "region", "group")
  )
  colnames(plot_defs) <- c("compartment", "treat_compare", "analysis")
  
  # 
  # print(comparisons)
  ## plots for all panels
  
  # scatter plots 
  # cannot use a for-loop...need to use Map or lapply
  lapply(1:nrow(plot_defs), function(i) {
    comp <- as.character(plot_defs[i, "compartment"])
    treat <- as.character(plot_defs[i, "treat_compare"])
    analysis <- as.character(plot_defs[i, "analysis"])
    
    pltname <- sprintf("panel%s_scatterplot_%s_%s", 
                       analysis, comp, treat)
    output[[pltname]] <- renderPlot({
      annot_df_name <- paste0(c("gene", analysis, "df"), collapse="_")
      annot_df <- annotation_dfs[[annot_df_name]]
      
      # if symbol or region plot gene symbols
      if (analysis == "group") {
        plot_scatter(comp, treat, annot_df, annotate_genes = FALSE)
      } else {
        plot_scatter(comp, treat, annot_df)
      }
      })
  })
  
  ## plots for select subsets
  # symbol only 
  # boxplots across samples
  lapply(unique(plot_defs$compartment), 
         function(comp) {
         
           comp <- as.character(comp)
           
           pltname <- sprintf("panelsymbol_boxplot_static_%s_all", comp)
           
           output[[pltname]] <- renderPlot({
             plot_boxplot(comp, annotation_dfs$gene_symbol_df)
           })
           
         })
  
    
  # group only
  lapply(unique(plot_defs$compartment),
         function(comp) {
           comp <- as.character(comp)
           
           pltname <- sprintf("panelgroup_violinplot_static_%s_all", comp)
           
           output[[pltname]] <- renderPlot({
             plot_violin_log2FC(comp, annotation_dfs$gene_group_df)
           })
           
         })

  
  
  #output$panelsymbol_lysoIP_plot_scatter_L.nigericin_L <- renderPlot({
  #  plot_interactive_scatter("lysoIP", "L.nigericin_L", 
  #                           annotation_dfs$gene_symbol_df)
  #})
  
  #output$panelsymbol_lysoIP_plot_scatter_L.CL097_L <- renderPlot({
  #  plot_interactive_scatter("lysoIP", "L.CL097_L", 
  #                           annotation_dfs$gene_symbol_df)
  #})
  
  #output[["panelsymbol_lysoIP_plot_scatter_L.CL097_L.nigericin"]] <- renderPlot({
  #  plot_interactive_scatter("lysoIP", "L.CL097_L.nigericin", 
  #                           annotation_dfs$gene_symbol_df)
  #})
  
  
  # Heatmap outputs
  # region
  lapply(unique(plot_defs$compartment),
         function(comp) {
           comp <- as.character(comp)
           
           pltname <- sprintf("panelregion_heatmap_static_%s_all", comp)
           print(pltname)
           print("trying to render heatmap")
           output[[pltname]] <- renderPlot({
             plot_heatmap(comp, annotation_dfs$gene_region_df)
           })
         })
  
  # group
  lapply(unique(plot_defs$compartment),
         function(comp) {
           comp <- as.character(comp)
           
           pltname <- sprintf("panelgroup_heatmap_static_%s_all", comp)
           output[[pltname]] <- renderPlot({
             plot_heatmap(comp, annotation_dfs$gene_group_df)
           })
         })
  
  
#  output$heatmap_region <- renderPlot({plot_heatmap(annotation_dfs$gene_region_df, "lysoIP")})
  
#  output$heatmap_group <- renderPlot({plot_heatmap(annotation_dfs$gene_group_df, "lysoIP")})
  
  
  # Table outputs
  lapply(c("symbol", "region", "group"),
         function(analysis) {
           df_name <- sprintf("gene_%s_df", analysis)
           tbl_name <- sprintf("table_%s", analysis)
           
           
           msmt_cols <- c(sprintf("log2FC_%s", unique(processed_data[['comparisons']]$comparison)),
                          sprintf("qval_%s", unique(processed_data[['comparisons']]$comparison)),
                          sprintf("%s_Intensity_Mean", unique(processed_data[['design']]$condition))        
                                 )
           
           output[[tbl_name]] <- DT::renderDataTable({DT::datatable(get_wide_data(annotation_dfs[[df_name]], msmt_cols),
                                                                    options=list(scrollX=T))
           })
           
         })
  
  print("output options")
  print(outputOptions(output))
  
}

# Run the application 
shinyApp(ui = ui, server = server)

