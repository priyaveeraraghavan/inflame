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

project_dir <- "/Users/priyaveeraraghavan/Documents/NLRP3_Bobby/Proteomics_Visualization"
load(file.path(project_dir, "data/processed_proteomics_data.RData"))
metadata_file <- file.path(project_dir, "Proteomics_Visualization/data/designNoCTR_msSTATS.csv")
input_data_file <- file.path(project_dir, "Proteomics_Visualization/data/all_stats_results_Lyso_1peptide.csv")

accent_colors <- c("red", "blue")

load_metadata <- function(meta_file) {
  meta_df <- read.csv(meta_file) %>%
    subset(condition != "Empty" | !is.na(experiment_type)) %>%
    rowwise() %>%
    mutate(sample_name=sprintf("%s_%s_%s", experiment_type, condition, replicate))
  
  return(meta_df)
}
meta_df <- load_metadata(metadata_file)

load_data <- function(input_file) {
  
  input_df <- read.csv(input_file)
  
  # filter out the "#" 
  input_df <- input_df %>% 
    subset(!grepl("#", reference) & nchar(gene_symbol) > 0) 
  
  # add negative log10 q values
  for (qval in colnames(input_df)[grepl("q.val", colnames(input_df))]) {
    input_df[,sprintf("negLog10%s", qval)] <- -log10(input_df[,qval])
  }
  
  # just for now with only lyso data
  gene_symbols <- input_df$gene_symbol
  references <- input_df$reference
  input_df$gene_symbol <- NULL
  input_df$reference <- NULL
  colnames(input_df) <- sprintf("lysoIP_%s", colnames(input_df))
  input_df$gene_symbol <- gene_symbols
  input_df$reference <- references
  return(input_df)
  
}

input_df <- load_data(input_data_file)

# get a list of comparisons based on the metadata and on what measurements
# exist in the input data frame
# searches for log2FC columns
get_comparisons <- function(meta_df, input_df) {
  poss_comparisons <- data.frame(
    expand.grid(list(condition1=unique(meta_df$condition), 
                     condition2=unique(meta_df$condition),
                     experiment_type=unique(meta_df$experiment_type))),
    stringsAsFactors = FALSE) %>%
    rowwise() %>% 
    mutate(comparison_name=sprintf("%s_%s", 
                                   as.character(condition1),
                                   as.character(condition2))) %>%
    rowwise() %>%
    mutate(log2FC_comparison_name=sprintf("%s_log2FC_%s", 
                                          as.character(experiment_type),
                                          comparison_name))
  
  # subset for measured columns
  valid_comparisons <- poss_comparisons %>%
    subset(log2FC_comparison_name %in% colnames(input_df))
  
  return(valid_comparisons)
  
}

comparison_df <- get_comparisons(meta_df, input_df)


## UI-side plot definitions
define_static_plot <- function(panel, plot_type, organelle, ip_comparison) {
  
  return(box(title=sprintf("%s %s", organelle, gsub("_", "/", ip_comparison)), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_%s_static_%s_%s", panel, plot_type, organelle, ip_comparison), 
                    width = 250, height=250)))
}

define_static_scatter <- function(panel, organelle, ip_comparison) {
  return(box(title=sprintf("%s %s", organelle, gsub("_", "/", ip_comparison)), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_scatterplot_%s_%s", panel, organelle, ip_comparison), 
                    width = 250, height=250)))
}

define_interactive_scatter <- function(panel, organelle, ip_comparison) {
  return(box(title=sprintf("%s %s", organelle, gsub("_", "/", ip_comparison)), width = NULL, height = 300,
  plotOutput(outputId = sprintf("%s_scatterplot_%s_%s", panel, organelle, ip_comparison), 
                    width = 250, height=250, 
                    brush = brushOpts(id=sprintf("%s_plot_brush", panel), resetOnNew = TRUE))))
}


## Helper plotting functions for Server side
plot_heatmap <- function(annot_df, experiment) {
  
  meta_subs <- subset(meta_df, experiment_type == experiment)
  
  heatmap_df <- annot_df[,meta_subs$sample_name]
  
  # normalize as prop of max
  #rowmaxs <- apply(heatmap_df, 1, max)
  heatmap_df <- data.frame(sweep(heatmap_df, MARGIN = 1, FUN = "/", STATS = rowMeans(heatmap_df)))
  
  # make sure rows and cols named
  rownames(heatmap_df) <- annot_df$reference
  colnames(heatmap_df) <- meta_subs$sample_name
  
  # define gene group for annotation
  annot_row <- data.frame(gene_group=annot_df$gene_group)
  rownames(annot_row) <- annot_df$reference
  
  # define sample group for annotation
  annot_col <- data.frame(meta_subs[,c("condition", "experiment_type")])
  rownames(annot_col) <- meta_df$sample_name
  pheatmap(heatmap_df, 
           annotation_row = annot_row,
           annotation_col = annot_col)
  
}

plot_scatter <- function(experiment, comparison, annot_df, annotate_genes=TRUE) {
  xcol <- paste0(c(experiment, "log2FC", comparison), collapse="_")
  ycol <- paste0(c(experiment, "negLog10q.val", comparison), collapse="_")
  
  print(xcol)
  print(ycol)
  
  plt <- ggplot(input_df, aes_string(xcol, ycol)) +
    geom_point(show.legend=FALSE)  +
    geom_point(data=annot_df, 
               mapping=aes_string(xcol, ycol, color="gene_group"),
               size=2,
               show.legend = FALSE) +
    theme_classic() +
    geom_vline(xintercept=0.5, linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_hline(yintercept=log2(0.05), linetype="dashed")
  
  if (annotate_genes) {
    plt <- plt + geom_label_repel(data=annot_df, 
                     mapping=aes_string(xcol, ycol, color="gene_group",
                                        label="gene_symbol"),
                     size=3,
                     show.legend = FALSE)
  }
  return(plt)
}

# for gene symbols
plot_boxplot <- function(experiment, annot_df) {
  
  sample_names <- subset(meta_df, experiment_type == experiment)$sample_name
  
  annot_df_subs <- annot_df[,c("gene_symbol", sample_names)] %>%
    pivot_longer(contains(experiment), 
                 names_to = "sample_name", values_to="intensity") %>%
    inner_join(meta_df, by="sample_name")
  
  # plot
  ggplot(annot_df_subs, aes(condition, intensity)) + 
    geom_boxplot() + 
    geom_point() +
    theme_classic()
  
}

# for gene group 
# violin plot log2FC
plot_violin_log2FC <- function(experiment, annot_df) {
  
  annot_df_subs <- annot_df %>% 
    dplyr::select((contains(experiment) & contains("log2FC")) | matches("gene_group")) %>%
    pivot_longer(contains("log2FC"), names_to="comparison", values_to="log2FC") %>%
    rowwise() %>%
    mutate(comparison=gsub("_", "/\n", gsub("_log2FC_", "", gsub(experiment, "", comparison))))
  
  # plot
  ggplot(annot_df_subs,
         aes(comparison, log2FC, color=gene_group)) + 
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
    theme_classic() +
    theme(legend.position="bottom") +
    ggtitle(experiment)
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
                                column(width=4,
                                       lapply(subset(comparison_df, experiment_type == "lysoIP")$comparison_name,
                                              function(comp) { define_static_scatter("panelsymbol", "lysoIP", comp)})
                                )
                              ),
                              fluidRow(
                                define_static_plot("panelsymbol", "boxplot", "lysoIP",  "all")
                                ), # in future will include the others
                              
                              fluidRow(
                                box(width = 12, DT::dataTableOutput("table_symbol")))
                      ),
    tabPanel("Gene Region Explorer",
             fluidRow(
               column(width=4,
                      lapply(subset(comparison_df, experiment_type == "lysoIP")$comparison_name,
                             function(comp) { define_interactive_scatter("panelregion", "lysoIP", comp)})
               )
             ),
             fluidRow(
               box("Gene Abundance Heatmap", width=12, plotOutput("heatmap_region"))
             ),
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
               column(width=4,
                      lapply(subset(comparison_df, experiment_type == "lysoIP")$comparison_name,
                             function(comp) { define_static_scatter("panelgroup", "lysoIP", comp)})
               )
             ),
             fluidRow(
               define_static_plot("panelgroup", "violinplot", "lysoIP",  "all")
             ), # in future will include the others
             fluidRow(
               box("Gene Abundance Heatmap", width=12, plotOutput("heatmap_group"))
             ),
             fluidRow(
               box(width = 12, DT::dataTableOutput("table_group"))
             )
    )
  )
)

                

server <- function(input, output) {
  
  empty_df <- input_df %>% 
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
      annotation_dfs$gene_symbol_df <- subset(input_df, gene_symbol %in% gene_symbols) %>%
        rowwise() %>% mutate(gene_group = "Gene_Selection")
    }
    
    print(annotation_dfs)
  })
  
  observeEvent(input$panelregion_plot_brush, ignoreNULL = TRUE, {
    brush_output_id <- input$panelregion_plot_brush$outputId
    brush_compare <- gsub("panelregion_scatterplot_", "", brush_output_id)
    experiment_type <- strsplit(brush_compare, "_")[[1]][1]
    comparison_name <- paste0(strsplit(brush_compare, "_")[[1]][2:3], collapse="_")
    
    brush_xcol <- paste0(c(experiment_type, "log2FC", comparison_name), collapse="_")
    brush_ycol <- paste0(c(experiment_type, "negLog10q.val", comparison_name), collapse="_")
    
    
    annotation_dfs$gene_region_df <- input_df %>% 
      dplyr::filter(.data[[brush_xcol]] < input$panelregion_plot_brush$xmax &
                      .data[[brush_xcol]] > input$panelregion_plot_brush$xmin &
                      .data[[brush_ycol]] < input$panelregion_plot_brush$ymax & 
                      .data[[brush_ycol]] > input$panelregion_plot_brush$ymin) %>%
      rowwise() %>% 
      mutate(gene_group="Regional_Selection")
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
    missing_ids <- setdiff(gene_gp_df$gene_symbol, input_df$gene_symbol)
    num_unique_ids <- length(unique(gene_gp_df$gene_symbol))
    map_pct <- (num_unique_ids-length(missing_ids))/length(num_unique_ids)*100
    
    # merge with input_df
    annotation_dfs$gene_group_df <- input_df %>% inner_join(gene_gp_df)
    
  })
  
  comparisons <- expand.grid(c("lysoIP"),
                             c("L.nigericin_L", "L.CL097_L", "L.CL097_L.nigericin"),
                             c("symbol", "region", "group"))
  colnames(comparisons) <- c("experiment", "treat_compare", "analysis")

  print(comparisons)
  ## plots for all panels
  
  # scatter plots 
  # cannot use a for-loop...need to use Map or lapply
  lapply(1:nrow(comparisons), function(i) {
    expt <- as.character(comparisons[i, "experiment"])
    treat <- as.character(comparisons[i, "treat_compare"])
    analysis <- as.character(comparisons[i, "analysis"])
    
    pltname <- sprintf("panel%s_scatterplot_%s_%s", 
                       analysis, expt, treat)
    print(pltname)
    output[[pltname]] <- renderPlot({
      annot_df_name <- paste0(c("gene", analysis, "df"), collapse="_")
      annot_df <- annotation_dfs[[annot_df_name]]
      
      # if symbol or region plot gene symbols
      if (expt == "group") {
        plot_scatter(expt, treat, annot_df, annotate_genes = FALSE)
      } else {
        plot_scatter(expt, treat, annot_df)
      }
      })
  })
  
  ## plots for select subsets
  # symbol only 
  # boxplots across samples
  lapply(unique(comparisons$experiment), 
         function(expt) {
         
           expt <- as.character(expt)
           
           pltname <- sprintf("panelsymbol_boxplot_static_%s_all", expt)
           
           output[[pltname]] <- renderPlot({
             plot_boxplot(expt, annotation_dfs$gene_symbol_df)
           })
           
         })
  
    
  # group only
  lapply(unique(comparisons$experiment), 
         function(expt) {
           expt <- as.character(expt)
           
           pltname <- sprintf("panelgroup_violinplot_static_%s_all", expt)
           
           output[[pltname]] <- renderPlot({
             print("plotting violin")
             plot_violin_log2FC(expt, annotation_dfs$gene_group_df)
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
  output$heatmap_region <- renderPlot({plot_heatmap(annotation_dfs$gene_region_df, "lysoIP")})
  
  output$heatmap_group <- renderPlot({plot_heatmap(annotation_dfs$gene_group_df, "lysoIP")})
  
  
  # Table outputs
  lapply(c("symbol", "region", "group"),
         function(analysis) {
           df_name <- sprintf("gene_%s_df", analysis)
           tbl_name <- sprintf("table_%s", analysis)
           
           output[[tbl_name]] <- DT::renderDataTable({DT::datatable(annotation_dfs[[df_name]],
                                                                    options=list(scrollX=T))
           })
           
         })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

