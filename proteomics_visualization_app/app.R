library(DT)
library(ggrepel)
library(ggpubr)
library(shiny)
library(tidyverse)
library(shinyjs)
library(shinydashboard)
library(RColorBrewer)
library(pheatmap)
library(bslib)

# Setup data
load("processed_proteomics_data.RData")

get_dataset_from_type <- function(dataset_type) {
  # quick way to get all original data for a type of data
  if (dataset_type == "ip_data") {
    return(fixed_processed_ip_data)
  } else if (dataset_type == "apex_sp_data") {
    return(fixed_processed_apex_singleplex_data)
  } else if (dataset_type == "apex_time_data") {
    return(fixed_processed_apex_timecourse_data) 
  } else {
    return(NULL)
  }
}


accent_colors <- c("#999999",  "#D55E00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2","#E69F00", "#F0E442")
get_color_intensity <- function(color_name, extreme) {
  fc <- colorRampPalette(c(color_name, extreme))(3)
  
  return(fc[2])  
}
accent_colors_shade <- sapply(accent_colors, function(col) get_color_intensity(col, "black"))
accent_colors_tint <- sapply(accent_colors, function(col) get_color_intensity(col, "white"))

##################
### Data Wrangling
##################

update_data_annotations <- function(dataset_info, 
                                    new_gene_annot_df=data.frame(reference=character(0),
                                                                 gene_group=character(0)), 
                                    filter_data=FALSE) {
  ## Typically used after a reactive event
  # should do it separately for each dataset type
  # new_gene_annot_df must have columns: gene_symbol, group_annotation
  # others will be dropped
  
  
  # get names of objs in processed data that will have gene names
  gene_objs <- names(dataset_info)[!grepl("design|comparisons|gene_annotations",names(dataset_info))]
  # merge, left joining to retain non-annotated genes
  if (nrow(new_gene_annot_df) > 0) {
    new_gene_annot_df <- new_gene_annot_df[,c("reference", "gene_group")]
    
    annot_dataset <- lapply(gene_objs, 
                            function(g) dataset_info[[g]] %>% 
                              left_join(new_gene_annot_df, by=c("reference")))
    
  } else {
    annot_dataset <- lapply(gene_objs, 
                            function(g) dataset_info[[g]] %>% rowwise() %>%
                              mutate(gene_group=NA)
    )
  }
  
  # decide if filter or not
  if (filter_data) {
    annot_dataset <- lapply(annot_dataset, 
                            function(g) {
                              g$show <- g$filtered
                              return(g)
                            })
  } else {
    annot_dataset <- lapply(annot_dataset, 
                            function(g) {
                              g$show <- TRUE
                              return(g)
                            })
  }
  
  # must rename to retain info about what is in each sub-dataset
  names(annot_dataset) <- gene_objs
  
  return(annot_dataset)
}

get_measurement_columns <- function(df, measurement_regexes) {
  df_cols <- colnames(df)
  msmt_cols <- unlist(sapply(measurement_regexes, 
                      function(rxp) df_cols[grepl(rxp, df_cols)]))
  return(unique(msmt_cols))
}

get_single_data_table_output <- function(df, measurement_regexes, 
                                         gene_identifiers=c("reference", "gene_symbol_species"), 
                                         reference_filter=NULL,
                                         filter_data=FALSE) {
  
  msmt_cols <- get_measurement_columns(df, measurement_regexes)

  df_select <- df[,c(msmt_cols, gene_identifiers, "filtered")]

  if (filter_data) {
    df_select <- df_select %>%
      # subset on filtered
      subset(filtered) 
  }
  
  # remove filtered
  df_select <- df_select %>%
    dplyr::select(-filtered)
  
  if (!is.null(reference_filter)) {
    df_select <- df_select %>%
      subset(reference %in% reference_filter)
  }
  
  return(df_select)
}

get_wide_multiple_data_table_output <- function(ds_obj_lst, ds_names, measurement_regexes, 
                                                gene_identifiers=c("reference", "gene_symbol_species"),
                                                reference_filter=NULL,
                                                filter_data=FALSE) {

  msmt_long <- do.call(rbind,
                       lapply(ds_names,
                              function(ds_name) {
                               
                                return(get_single_data_table_output(ds_obj_lst[[ds_name]], 
                                                             measurement_regexes,
                                                             reference_filter = reference_filter,
                                                             filter_data = filter_data) %>%
                                  rowwise() %>%
                                  mutate(dataset_name=ds_name))
                              }))
  
  # get cols to pivot
  column_names <- setdiff(colnames(msmt_long), c(gene_identifiers, "dataset_name"))
  print(intersect(colnames(msmt_long), gene_identifiers))
  msmt_wide <- msmt_long %>% pivot_wider(id_cols=c("reference", "gene_symbol_species"), # should work but doesn'ttidyr::all_of(gene_identifers), 
                                         names_from=dataset_name,
                                         values_from=column_names)
  
  return(msmt_wide)
}


get_dataset_measurements <- function(annotated_data, dataset_names, column_names) {
  
  # make sure that column_names also includes boolean of "show"
  column_names <- unique(c(column_names, "show"))
  
  # grab and bind the relevant columns
  # important that all dataframes have the same column_names!
  df <- do.call(rbind,
                lapply(dataset_names,
                       function(name) annotated_data[[name]][,column_names] %>%
                         rowwise() %>%
                         mutate(dataset_name=name)
                )
  )
  
  # only return valid rows from data
  # for example, with IP data this is a filtration bool w.r.t. 
  # background enrich of organelle itself
  df <- subset(df, show)
  
  return(df)
  
}

get_wide_data <- function(processed_data, dataset_names, column_names, filter_bool=FALSE) {
  # function that gets the relevant information for organelles and WCP for a gene set
  # displays in wide format
  
  long_df <- get_dataset_measurements(processed_data, dataset_names, c("reference", "gene_symbol", "gene_group", column_names))
  
  # remove nas in gene_group
  long_df <- long_df %>% subset(!is.na(gene_group))
  
  wide_df <- long_df %>% pivot_wider(id_cols=c(reference, gene_symbol, gene_group),
                                     names_from=dataset_name,
                                     values_from=column_names)
  
  
  return(wide_df)
  
}    

get_long_intensities_time_lm <- function(lm_res_df) {
  # get intensities in long format with design info
  intensity_cols <- colnames(lm_res_df)[grepl("Intensity|estimate|qvalue", colnames(lm_res_df))]
  df <- lm_res_df[,c("reference", intensity_cols)] %>%
    pivot_longer(!reference & !contains("estimate") & !contains("qvalue"), names_to="sample_name", values_to="true_intensity")
  
  # extract design info
  df$condition <- ifelse(grepl("NLRP", df$sample_name), "NLRP3", "P4C")
  df$timepoint <- sapply(df$sample_name, 
                         function(sname) as.numeric(gsub("CTR|NLRP", "", strsplit(sname, "_")[[1]][1])))
  df$predicted_intensity <- df$estimate_Intercept + 
    df$estimate_conditionNLRP*(df$condition == "NLRP3") + df$estimate_timepoint*df$timepoint +
    df$estimate_timepoint.conditionNLRP*(df$condition == "NLRP3")*df$timepoint
  df$is_sig_interaction <- df$qvalue_timepoint.conditionNLRP < 0.05
  
  return(df)
}

######################
### Plotting Functions
######################
cleanup_axes_names <- function(ax_name) {
  ax_name <- strsplit(gsub("vs", " vs ", ax_name), "_")[[1]]
  ax_name <- paste(ax_name, collapse="\n")
  ax_name
}

plot_heatmap <- function(annotated_data, ds_name, design_info, gene_group_color_map=NULL) {
  # compartment_name: string compartment 
  # gene_group_df: df of ['reference', 'gene_group']
  

  # column_name in design refers to the column in the true df 
  # containing either mean or replicate-level intensity
  design_subs <- subset(design_info, 
                        !grepl("Mean", column_name) & dataset_name == ds_name)
  
  # design_cols remove "intensity"
  design_col_rename <- gsub("_Intensity", "", design_subs$column_name)
  intensity_subset_df <- annotated_data[[ds_name]] %>% 
    subset(show & !is.na(gene_group)) 
  
  # if zero size then return empty plot
  if (nrow(intensity_subset_df) == 0) {
    return(ggplot())
  } else if (nrow(intensity_subset_df) <= 2) {
    cluster_rows_bool <- FALSE
  } else {
    cluster_rows_bool <- TRUE
  }
  
  intensity_subset_df <- 
    intensity_subset_df[,  c("reference", "gene_group", design_subs$column_name)] 
  
  # make unique
  intensity_subset_df <- unique(intensity_subset_df)
  heatmap_df <- intensity_subset_df[,  design_subs$column_name]
  
  # normalize as prop of max
  #rowmaxs <- apply(heatmap_df, 1, max)
  heatmap_df <- data.frame(sweep(heatmap_df, MARGIN = 1, FUN = "/", STATS = rowMeans(heatmap_df)))
  
  # make sure rows and cols named
  rownames(heatmap_df) <- intensity_subset_df$reference
  colnames(heatmap_df) <- design_col_rename
  
  # define gene group for annotation
  annot_row <- data.frame(gene_group=intensity_subset_df$gene_group) #%>% drop_na()
  rownames(annot_row) <- intensity_subset_df$reference
  
  # define sample group for annotation
  annot_col <- data.frame(condition=design_subs$condition)
  rownames(annot_col) <- design_col_rename
  
  # get colors
  # for the heatmap itself
  heatmap_color <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                                     "Greys")))(100)
  # for annotations
  if (is.null(gene_group_color_map)) {
    gene_groups <- unique(annot_row$gene_group)
    n_gene_groups <- length(gene_groups)
    gene_group_color_map <- accent_colors[2:(n_gene_groups+1)]
    names(gene_group_color_map) <- gene_groups
  } else {
    gene_group_color_map <- gene_group_color_map[2:length(gene_group_color_map)]
  }
  
  conditions <- unique(annot_col$condition)
  n_conditions <- length(conditions)
  condition_colors <- accent_colors_tint[2:(n_conditions+1)]
  names(condition_colors) <- conditions
  
  annot_colors <- list(gene_group=gene_group_color_map,
                       condition=condition_colors)
  
  pheatmap(heatmap_df, color=heatmap_color,
           annotation_colors=annot_colors,
           annotation_row = annot_row,
           annotation_col = annot_col,
           cluster_rows = cluster_rows_bool)
  
}

plot_scatter <- function(annotated_data, ds_name, design_info, xcol, ycol, 
                         annotate_genes=TRUE, max_genes_annot=10,
                         gene_group_color_map=NULL) {

  data_full <- data.frame(annotated_data[[ds_name]] %>% 
                            subset(show))
  

  data_subset <- data.frame(data_full %>%
                              subset(!is.na(gene_group)))
  
  
  plt <- ggplot(data_full, aes_string(xcol, ycol),color='black') + 
    geom_point(alpha=0.5) 
  
  if (nrow(data_subset) > 0) {
    plt <- plt +
      geom_point(data=data_subset, 
                 mapping=aes_string(xcol, ycol, color="gene_group"),
                 size=2,
                 show.legend = FALSE)
  }
  
  plt <- plt +
    theme_classic() +
    geom_vline(xintercept=0.5, linetype="dashed") +
    geom_vline(xintercept=-0.5, linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") + 
    xlab(cleanup_axes_names(xcol)) + 
    ylab(cleanup_axes_names(ycol)) 
  
  if (annotate_genes & nrow(data_subset) > 0) {
    data_subset$gene_score <- abs(data_subset[xcol])*data_subset[ycol]
    to_annotate <- data_subset %>%
      arrange(desc(gene_score)) %>%
      head(max_genes_annot)
    
    plt <- plt + geom_label_repel(data=to_annotate, 
                                  mapping=aes_string(xcol, ycol, color="gene_group",
                                                     label="gene_symbol_species"),
                                  size=3,
                                  show.legend = FALSE)
  }
  
  if (!is.null(gene_group_color_map)) {
    plt <- plt + scale_color_manual(values=gene_group_color_map)
  }
  
  
  return(plt)
  
}

# plot boxplot of per-replicate intensities
plot_boxplot_intensity <- function(annotated_data, ds_name, design_info, gene_group_color_map=NULL) {
  design_subset <- subset(design_info,
                          dataset_name == ds_name & replicate != "Mean")
  
  column_names <- design_subset$column_name
  
  data_subset <- data.frame(annotated_data[[ds_name]] %>% 
                              subset(show & !is.na(gene_group)))
  
  if (nrow(data_subset) > 0) {
    
    annot_df_subs <- data_subset[,c("reference", column_names)] %>%
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
plot_violin_msmt <- function(annotated_data, ds_name, 
                             msmt_col, gene_group_color_map=NULL) { 
  
  data_all <- get_dataset_measurements(annotated_data, ds_name, 
                                          unique(c("reference", "gene_group", msmt_col))) 
  
  data_subset <- data_all %>%
    subset(!is.na(gene_group))
  

  
  if (length(unique(data_subset$gene_group)) == 1) {
    data_all$gene_group <- "background"
    data_subset <- rbind(data_subset, data_all)
  } 
  
  # use gene group sizes to determine whether violin or boxplot
  gene_group_sizes <- (data_subset %>%
                         group_by(gene_group) %>%
                         summarize(number=n()))$number
  if (max(gene_group_sizes) < 30) {
    plot_type <- "box"
  } else {
    plot_type <- "violin"
  }
  
  annot_df_subs <- data_subset[,c(msmt_col, "gene_group")] %>%
    pivot_longer(!gene_group, names_to="comparison", values_to="measurement") %>%
    unique()
  
  # modify the x axis names so they fit
  annot_df_subs <- annot_df_subs %>%
    rowwise() %>%
    mutate(comparison=gsub("vs", "\nvs\n",
                           gsub("_", "", 
                                gsub("log2FC", "", comparison))))
  
  
  # plot
  plt <- ggplot(annot_df_subs)
  if (plot_type == "violin") {
  plt <-  ggplot(annot_df_subs,
                aes(comparison, measurement, color=gene_group)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
    theme_classic() +
    theme(legend.position="bottom") +
    ggtitle(ds_name)
  } else {
      plt <- ggplot(annot_df_subs, aes(comparison, measurement, color=gene_group)) +
        geom_boxplot() +
        geom_point(position=position_jitterdodge(jitter.width=0.1)) +
        theme_classic() +
        theme(legend.position="bottom") +
        ggtitle(ds_name)
  }
  
  if (!is.null(gene_group_color_map)) {
    plt <- plt + scale_color_manual(values=gene_group_color_map)
  }
  
  return(plt)
}

plot_timecourse_linear_model <- function(lm_res_df,
                                         gene_group_color_map=NULL) {
  # grab relevant/selected data
  selected_data <- lm_res_df %>% subset(!is.na(gene_group) & show)
  
  
  if (nrow(selected_data) == 0) {
    return(ggplot())
  }
  # if there are more than 4 genes, then plot a violin instead
  if (length(unique(selected_data$gene_symbol_species)) == 1) {
    # put in long format for plotting
    long_intensities <- get_long_intensities_time_lm(selected_data)
    # plot
    color_map <- accent_colors[1:2]
    names(color_map) <- c(FALSE, TRUE)
    plt <- long_intensities %>%
      ggplot() +
      geom_point(mapping=aes(timepoint, true_intensity, color=is_sig_interaction, shape=condition), size=1, position=position_jitterdodge(jitter.width=0.5, dodge=2)) + 
      geom_line(mapping=aes(timepoint, predicted_intensity, color=is_sig_interaction, linetype=condition), alpha=0.5, linewidth=1) +
      theme_classic() +
      scale_color_manual(values=color_map) + 
      facet_wrap(vars(reference),nrow=2, scales="free_y") +
      ylab("Log2 Normalized MS Intensity") +
      xlab("Time (min)") +
      theme(legend.position = "bottom")
  } else {
    
    # rename columns
    lm_res_df <- lm_res_df %>% dplyr::rename(b_NLRP3=estimate_conditionNLRP, 
                                             b_time=estimate_timepoint,
                                             b_time_NLRP3=estimate_timepoint.conditionNLRP)
    
    # pivot coeffs longer
    long_coefs_all <- lm_res_df %>% subset(show) %>%
      pivot_longer(contains("b_") & !contains("Intercept"), 
                   names_to="coefficient_name", values_to="estimate") %>%
      rowwise() %>% mutate(gene_group="background")
    
    
    long_coefs_selected <- lm_res_df %>% subset(!is.na(gene_group) & show) %>%
      pivot_longer(contains("b_") & !contains("Intercept"), 
                   names_to="coefficient_name", values_to="estimate") 
    
    long_coefs <- rbind(long_coefs_selected, long_coefs_all)
    
    gene_group_levels <- c("background", 
                           unique(long_coefs_selected$gene_group))
    if (length(gene_group_levels) > length(accent_colors)) {
      print("too many colors")
    }
    
    long_coefs$gene_group <- factor(long_coefs$gene_group, 
                                    levels=gene_group_levels)
    
    plt <- ggplot(long_coefs)
    if (length(unique(selected_data$reference)) <= 20) {
      plt <- plt + geom_boxplot(mapping=aes(gene_group, estimate, color=gene_group)) +
        geom_point(data=long_coefs_selected, 
                   mapping=aes(gene_group, estimate, color=gene_group), 
                   position=position_jitter(width=0.2, height = 0)) 
    } else {
      plt <- plt + geom_violin(mapping=aes(gene_group, estimate, color=gene_group),
                               draw_quantiles = c(0.25, 0.5, 0.75)) 
      
    }
    if (!is.null(gene_group_color_map)) {
      gene_group_color_map <- accent_colors[1:length(gene_group_levels)]
    }
    plt <- plt + 
      theme_classic() +
      theme(legend.position="None") +
      scale_color_manual(values=gene_group_color_map) +
      geom_hline(yintercept=0, color='black', linetype='dashed') +
      facet_wrap(vars(coefficient_name), ncol=1, scales="free") + 
      xlab("") +
      stat_compare_means(aes(x = gene_group, y = estimate), label = "p.signif", method = "anova", label.x.npc = "center") +
      coord_flip()
    
  }
  
  plt
  
  
}

#############
########## UI 
#############

## UI-side plot definitions
format_label <- function(label) {
  # replace underscores
  label <- gsub("_", " ", label)
  
  # capitalize
  label <- str_to_title(label)
  
  # vs gets spaces
  label <- gsub("vs", " vs ", label)
  
  # replace dots with : for interaction
  label <- gsub("\\.", ":", label)
  
  # fix cellLine
  label <- gsub("Cellline|cellline", "CellLine", label)
  
  # fix treatment names
  label <- gsub("cl097|Cl097", "CL097", label)
  label <- gsub("nigericin", "Nigericin", label)
  label <- gsub("control", "Control", label)
  label <- gsub("P4c|p4c", "P4C", label)
  
  # remove "any" or "all" at the end
  label <- gsub("Any|any|All$|all$", "", label)
  
  # change treat
  label <- gsub(" Treat ", " All Treatment ", label)
  
  return(label)
  
}
format_title <- function(labels, fontsize=15) {
  
  text_cleaned <- paste(sapply(labels, format_label), collapse="\n")
  
  return(h4(text_cleaned, style=sprintf('font-size:%spx;color:black;', fontsize)))
}

define_static_plot <- function(panel, plot_type, dataset_name, treatvsctrl) {
  
  return(box(title=format_title(c(dataset_name, treatvsctrl)), width = NULL, height = 300,
             plotOutput(outputId = sprintf("%s-%s_static-%s-%s", panel, plot_type, dataset_name, treatvsctrl),
                        height="250px")))
}

define_interactive_scatter <- function(panel, dataset_name, treatvsctrl) {
  return(box(title=format_title(c(dataset_name, treatvsctrl)), width = NULL, height = 300,
             plotOutput(outputId = sprintf("%s-scatterplot_reactive-%s-%s", panel, dataset_name, treatvsctrl), 
                        height="250px",
                        brush = brushOpts(id=sprintf("%s_plot_brush", panel), resetOnNew = TRUE))))
}


get_apex_sp_choices <- function() {
  gsub("vs", " vs ", 
       get_dataset_from_type("apex_sp_data")$comparisons$comparison)
}

ip_datasets <- unique(get_dataset_from_type("ip_data")$design$dataset_name)
ip_comparisons_df <- get_dataset_from_type("ip_data")$comparisons

sp_datasets <- c("all_bait_treat_single_plex", "cellLine.Nigericin_linear_model", 
                 "cellLine.CL097_linear_model", "CL097.Nigericin_linear_model")


timecourse_datasets <- c("apex_time_dataset_select", "nigericin_timecourse_linear_model", 
                         "cl097_timecourse_linear_model")

timecourse_heatmap_datasets <- c("nlrp3_nigericin_timecourse", "nlrp3_cl097_timecourse",
                                 "p4c_nigericin_timecourse", "p4c_cl097_timecourse")

## Define UI for application 
ui <- fluidPage(
  titlePanel("INFLAME"),
  tabsetPanel(id = "tabs",
    tabPanel("Organelle IP",
             fluidRow(
               box(width=5, selectInput("ip_input_type", "Choose Input Type",
                           choices = c("Single Gene Viewer", "Gene Region Explorer", "Gene Group Viewer"),
                           selected = "Single Gene Viewer")
               )
               ),
               fluidRow(
                 box(width=5,
                 selectInput("ip_filter", "Filter IP Dataset",
                             choices = c("Yes", "No"),
                             selected = "Yes")
                 )
               ),
               conditionalPanel(
                 condition = "input.ip_input_type == 'Single Gene Viewer'",
                 fluidRow(box(
                   width = 5,
                   textInput(
                     inputId = "ip__symbol",
                     label = "Gene Symbol",
                     value = "",
                     placeholder = "e.g., ACTB, E9Q5F4, ENSG00000075624"
                   )     # might also want to try textAreaInput for longer
                 )),
                 fluidRow(lapply(ip_datasets,
                                 function(ds_name)
                                   column(
                                     width = 2,
                                     lapply(subset(ip_comparisons_df, dataset_name == ds_name)$comparison,
                                            function(comparison) {
                                              define_static_plot("ip__panelsymbol", "scatterplot", ds_name, comparison)
                                            }),
                                   ))),
                 fluidRow(lapply(ip_datasets,
                                 function(ds_name)
                                   column(
                                     width = 2,
                                     define_static_plot("ip__panelsymbol", "boxplot", ds_name,  "all")
                                   ))),
                 fluidRow(box(width = 12, DT::dataTableOutput("ip__panelsymbol_table")))
               ),
               conditionalPanel(
                 condition = "input.ip_input_type == 'Gene Region Explorer'",
                 fluidRow(box(width=12, p("Select and drag to capture genes in a region of any scatterplot."))),
                 fluidRow(lapply(ip_datasets,
                                 function(ds_name)
                                   column(
                                     width = 2,
                                     lapply(subset(ip_comparisons_df, dataset_name == ds_name)$comparison,
                                            function(comparison) {
                                              define_interactive_scatter("ip__panelregion", ds_name, comparison)
                                            })
                                   ))),
                 lapply(ip_datasets,
                        function(ds_name)
                          fluidRow(box(
                            format_title("Gene Abundance Heatmap", ds_name),
                            width =
                              10,
                            plotOutput(sprintf(
                              "ip__panelregion_heatmap_static_%s_all", ds_name
                            ))
                          ))),
                 fluidRow(box(width = 12, DT::dataTableOutput("ip__panelregion_table")))
               ),
               conditionalPanel(
                 condition = "input.ip_input_type == 'Gene Group Viewer'",
                 fluidRow(
                   box(width=12, p("Please upload a .csv with columns [gene_id, gene_group]. Valid gene_ids are gene name, UNIPROT identifier, and ENSEMBL identifier."))
                 ),
                 fluidRow(
                   box(width = 5,
                       fileInput(
                     "ip__gene_group_input",
                     "Gene Groupings",
                     accept = ".csv",
                     buttonLabel = "Upload"
                       )
                   )
                 ),
                 fluidRow(lapply(ip_datasets,
                                 function(ds_name)
                                   column(
                                     width = 2,
                                     lapply(subset(ip_comparisons_df, dataset_name == ds_name)$comparison,
                                            function(comparison) {
                                              define_static_plot("ip__panelgroup", "scatterplot", ds_name, comparison)
                                            })
                                   ))),
                 fluidRow(lapply(ip_datasets,
                                 function(ds_name)
                                   column(
                                     width = 2,
                                     define_static_plot("ip__panelgroup", "violinplot", ds_name,  "all")
                                   ))),
                 lapply(ip_datasets,
                        function(ds_name)
                          fluidRow(box(
                            format_title("Gene Abundance Heatmap", ds_name),
                            width =
                              10,
                            plotOutput(sprintf(
                              "ip__panelgroup_heatmap_static_%s_all", ds_name
                            ))
                          ))),
                 fluidRow(box(width = 12, DT::dataTableOutput("ip__panelgroup_table")))
                 )
 
    ),
    tabPanel("APEX", id="apex",
             fluidRow(
               box(width = 5, selectInput("apex_input_type", "Choose Analysis Type",
                                        choices =  c("Single Gene Viewer", "Gene Region Explorer", "Gene Group Viewer"),
                                        selected = "Single Gene Viewer"))
             ),
             fluidRow(
               box(width = 5, selectInput("apex_experiment_type", "Choose Experiment",
                                        choices = c("Compare Conditions", "Timecourse"),
                                        selected = "Compare Conditions"))
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Compare Conditions'",
               fluidRow(box(width = 5, 
                            selectInput(inputId = "apex_sp_comparison", 
                                        label = "Choose APEX Comparison",
                                        choices=get_apex_sp_choices())))
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Compare Conditions' && input.apex_input_type == 'Single Gene Viewer'",
               fluidRow(
                 box(width = 5, textInput(inputId = "apex_sp__symbol", label = "Gene Symbol", 
                                        value = "", placeholder = "e.g., Nlrp3, Q8R4B8, ENSMUSG00000032691")     # might also want to try textAreaInput for longer
                 )),
               fluidRow(lapply(sp_datasets, 
                               function(ds_name) {
                                 column(width=3, 
                                        define_static_plot("apex_sp__panelsymbol",
                                                           "scatterplot",
                                                           ds_name,
                                                           "any"))
                               })
               ),
               fluidRow(box(width = 10, define_static_plot("apex_sp__panelsymbol", "boxplot", "all_bait_treat_single_plex", "all"))),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_sp__panelsymbol_table")))
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Compare Conditions' && input.apex_input_type == 'Gene Region Explorer'",
               fluidRow(box(width=12, p("Select and drag to capture genes in a region of any scatterplot."))),
               fluidRow(lapply(sp_datasets, 
                               function(ds_name) {
                                 # define the interactive scatter in order to select genes
                                 column(width=3, 
                                        define_interactive_scatter("apex_sp__panelregion", ds_name, "any")
                                 )
                                 
                               })
               ),
               fluidRow(lapply(sp_datasets, 
                               function(ds_name) {
                                 column(width=3, 
                                        define_static_plot("apex_sp__panelregion",
                                                           "violinplot",
                                                           ds_name,
                                                           "all"))
                               })
               ),
               fluidRow(box(format_title(c("Gene Abundance Heatmap")), 
                            width=10,
                            plotOutput(sprintf("apex_sp__panelregion_heatmap_static_all_bait_treat_single_plex_all"))
               )),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_sp__panelregion_table"))
               )
             ),
             conditionalPanel(                  
               condition = "input.apex_experiment_type == 'Compare Conditions' && input.apex_input_type == 'Gene Group Viewer'",
               fluidRow(box(width=12,
                            p("Please upload a .csv with columns [gene_id, gene_group]. Valid gene_ids are gene name, UNIPROT identifier, and ENSEMBL identifier."))),
               fluidRow(
                 box(width = 5,
                     fileInput(inputId = "apex_sp__gene_group_input", label = "Gene Groupings",
                           accept=".csv",
                           buttonLabel = "Upload"))
               ),
               fluidRow(lapply(sp_datasets, 
                               function(ds_name) {
                                 column(width=3, 
                                        define_static_plot("apex_sp__panelgroup",
                                                           "scatterplot",
                                                           ds_name,
                                                           "any"))
                               })
               ),
               fluidRow(lapply(sp_datasets, 
                               function(ds_name) {
                                 column(width=3, 
                                        define_static_plot("apex_sp__panelgroup",
                                                           "violinplot",
                                                           ds_name,
                                                           "all"))
                               })
                        
               ),
               fluidRow(box(format_title(c("Gene Abundance Heatmap")),
                            width=10,
                            plotOutput(sprintf("apex_sp__panelgroup_heatmap_static_all_bait_treat_single_plex_all")))
               ),
               #fluidRow(downloadButton('downloadData',"Download the data")),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_sp__panelgroup_table")))
             ),
             conditionalPanel(
               condition="input.apex_experiment_type == 'Timecourse'",
               fluidRow(box(width = 5, selectInput(inputId = "apex_time_dataset_select", 
                                                 label = "Choose APEX Timecourse Dataset",
                                                 choices=unique(get_dataset_from_type("apex_time_data")$design$dataset_name))),
                        box(width = 5, selectInput(inputId = "apex_time_comparison", 
                                                 label = "Choose APEX Timepoint Comparison",
                                                 choices=NULL)))
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Timecourse' && input.apex_input_type == 'Single Gene Viewer'",
               fluidRow(
                 box(width = 5, textInput(inputId = "apex_time__symbol", label = "Gene Symbol", 
                                        value = "", placeholder = "e.g., Nlrp3, Q8R4B8, ENSMUSG00000032691")     # might also want to try textAreaInput for longer
                 )),
               fluidRow(
                 lapply(timecourse_datasets, 
                        function(ds_name) {
                          column(width=4, define_static_plot("apex_time__panelsymbol",
                                                             "scatterplot",
                                                             ds_name,
                                                             "any"))
                        })
               ),
               #fluidRow(downloadButton('downloadData',"Download the data")),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_time__panelsymbol_table"))
               )
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Timecourse' && input.apex_input_type == 'Gene Region Explorer'",
               fluidRow(box(width=12, p("Select and drag to capture genes in a region of any scatterplot."))),
               fluidRow(
                 lapply(timecourse_datasets,
                        function(ds_name) {
                          column(width=4, define_interactive_scatter("apex_time__panelregion",
                                                                     ds_name,
                                                                     "any"))
                        })
               ),
               lapply(timecourse_heatmap_datasets,
                      function(ds_name) { fluidRow(box(format_title(c("Gene Abundance Heatmap", ds_name)), 
                                                       width=10, 
                                                       plotOutput(sprintf("apex_time__panelregion_heatmap_static_%s_all", ds_name))))
                      }
               ),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_time__panelregion_table"))
               )
             ),
             conditionalPanel(
               condition = "input.apex_experiment_type == 'Timecourse' && input.apex_input_type == 'Gene Group Viewer'",
               fluidRow(box(width=12,
                            p("Please upload a .csv with columns [gene_id, gene_group]. Valid gene_ids are gene name, UNIPROT identifier, and ENSEMBL identifier."))),
               fluidRow(
                 box(width = 5, fileInput(inputId = "apex_time__gene_group_input", label = "Gene Groupings",
                           accept=".csv",
                           buttonLabel = "Upload"))
               ),
               fluidRow(
                 lapply(
                   timecourse_datasets,
                   function(ds_name) {
                     column(width=4, define_static_plot("apex_time__panelgroup",
                                                        "scatterplot",
                                                        ds_name,
                                                        "any"))
                   })
               ),
               lapply(timecourse_heatmap_datasets,
                      function(ds_name) fluidRow(box(format_title(c("Gene Abundance Heatmap", ds_name)), 
                                                     width=10, 
                                                     plotOutput(sprintf("apex_time__panelgroup_heatmap_static_%s_all", ds_name))))
               ),
               fluidRow(
                 box(width = 12, DT::dataTableOutput("apex_time__panelgroup_table"))
               )
             )
    ),
    fluidRow(box(width = 5, downloadButton('downloadData',"Download the data")))
  )
)




######################
############### Server 
######################



server <- function(input, output, session) {
  

  
  #######################
  ## Reactives and Helper
  
  ## Reactive Values
  annotation_dfs <- reactiveValues(  
    # keeps track of the current gene or group of genes selected
    current_gene_group_df = data.frame(gene_group=character(0), 
                                       reference=character(0)),
    # color map for gene group so all plots are consistent
    gene_group_color_map = NULL,
    # updated ip data
    ip_data = update_data_annotations(get_dataset_from_type("ip_data"), filter_data=FALSE),
    # updated apex timecourse data
    apex_time_data = update_data_annotations(get_dataset_from_type("apex_time_data"), filter_data=FALSE),
    # updated apex data single time
    apex_sp_data = update_data_annotations(get_dataset_from_type("apex_sp_data"), filter_data=FALSE),
    # need to be able to download at any time
    data_table_output = data.frame()
  )
  
  ## Reactive Expressions

  get_current_panel <- reactive({
    
    input_type_choices <- c("symbol", "region", "group")
    names(input_type_choices) <- c("Single Gene Viewer", "Gene Region Explorer", "Gene Group Viewer")
    print("entered reactive")
    print(input$tabs)
    if (input$tabs == "Organelle IP") {
      gene_input <- input_type_choices[[input$ip_input_type]]
      return(paste0("ip__panel", gene_input))
      
      
    } else if (input$tabs == "APEX") {
      
      gene_input <- input_type_choices[[input$apex_input_type]]
      
      if (input$apex_experiment_type == "Timecourse") {
        return(paste0("apex_time__panel", gene_input))
        
      } else if (input$apex_experiment_type == "Compare Conditions") {
        return(paste0("apex_sp__panel", gene_input))
      }
    } else {
      return("")
    }
  })
  
  get_apex_time_comparison_choices <- reactive({
    possible_comparisons <- subset(
      get_dataset_from_type("apex_time_data")$comparisons,
      dataset_name == input$apex_time_dataset_select)$comparison
    
    return(possible_comparisons)
  })
  
  ip_data_table_output <- reactive({
    filter_status <- input$ip_filter == "Yes"
    
    ip_data <- get_dataset_from_type("ip_data")
    
    gene_ids <- c("reference", "gene_symbol_species")
    
    ds_names <- names(ip_data)[!grepl("comparisons|design", names(ip_data))]
    
    msmt_wide <- get_wide_multiple_data_table_output(ip_data, ds_names, 
                                        c("log2FC", "^qval", "Intensity_Mean"), 
                                        gene_identifiers=gene_ids,
                                        reference_filter=annotation_dfs$current_gene_group_df$reference,
                                        filter_data=filter_status) 
    
    # join to get gene groups
    msmt_wide <- msmt_wide %>%
      inner_join(annotation_dfs$current_gene_group_df, by=c("reference"))
    
    # order to display nicely
    gene_id_cols <- c(gene_ids, "gene_group")
    other_cols <- setdiff(colnames(msmt_wide), gene_id_cols)
    msmt_wide <- msmt_wide[,c(gene_id_cols, other_cols)]
    
    return(msmt_wide)

  })
  
  apex_time_data_table_output <- reactive({
    
    apex_time_data <- get_dataset_from_type("apex_time_data")
    
    gene_ids <- c("reference", "gene_symbol_species")
    
    # get wide dataset measurements, grouping by datasets that have identical column names
    wide_nigericin <- get_wide_multiple_data_table_output(
      apex_time_data, 
      c("nlrp3_nigericin_timecourse", "p4c_nigericin_timecourse"),
      c("log2FC", "^qval", "Intensity_Mean"), 
      gene_identifiers=gene_ids,
      reference_filter = annotation_dfs$current_gene_group_df$reference,
      filter_data=FALSE) 
    
    wide_cl097 <- get_wide_multiple_data_table_output(
      apex_time_data, 
      c("nlrp3_cl097_timecourse", "p4c_cl097_timecourse"),
      c("log2FC", "^qval", "Intensity_Mean"), 
      gene_identifiers=gene_ids,
      reference_filter = annotation_dfs$current_gene_group_df$reference,
      filter_data=FALSE) 
    
    wide_lm <- get_wide_multiple_data_table_output(
      apex_time_data,
      c("nigericin_timecourse_linear_model", "cl097_timecourse_linear_model"),
      c("estimate_timepoint", "estimate_conditionNLRP", "^qvalue_timepoint", "^qvalue_conditionNLRP"),
      gene_identifiers=gene_ids,
      reference_filter = annotation_dfs$current_gene_group_df$reference,
      filter_data=FALSE
    )
    
    # remove "timecourse" and "estimate" from colnames
    colnames(wide_nigericin) <- gsub("_timecourse", "", colnames(wide_nigericin))
    colnames(wide_cl097) <- gsub("_timecourse", "", colnames(wide_cl097))
    colnames(wide_lm) <- gsub("estimate", "b", colnames(wide_lm))
    
    # merge all different dsets
    msmt_wide <- purrr::reduce(list(wide_nigericin, wide_cl097, wide_lm),
                             function(x, y) {x %>% full_join(y, by=gene_ids)})
    
    # filter for current genes
    msmt_wide <- msmt_wide %>%
      inner_join(annotation_dfs$current_gene_group_df, by=c("reference"))
    
    # order to display nicely
    gene_id_cols <- c(gene_ids, "gene_group")
    other_cols <- setdiff(colnames(msmt_wide), gene_id_cols)
    msmt_wide <- msmt_wide[,c(gene_id_cols, other_cols)]
    
    return(msmt_wide)
    
  })
  
  apex_sp_data_table_output <- reactive({
    
    apex_sp_data <- get_dataset_from_type("apex_sp_data")
    
    gene_ids <- c("reference", "gene_symbol_species")
    
    print("getting interaction models")
    interaction_models <- lapply(
      c("cellLine.CL097", "cellLine.Nigericin", "CL097.Nigericin"),
      function(coef) {
        ds_name <- paste0(coef, "_linear_model")
        df <- apex_sp_data[[ds_name]]
        col_regex <- paste0(c("^qval_", "b_"), coef)
        return(get_single_data_table_output(df, col_regex,
                                     gene_identifiers = gene_ids,
                                     reference_filter = annotation_dfs$current_gene_group_df$reference,
                                     filter_data = FALSE))
      }
    )

    interaction_models_wide <- 
      purrr::reduce(interaction_models, 
                    function(x, y) {x %>% full_join(y, by=gene_ids)})
    
    print("getting single data table output")
    singleplex_wide <- get_single_data_table_output(
      apex_sp_data$all_bait_treat_single_plex,
      c("log2FC", "^qval", "Intensity_Mean"),
      gene_identifiers = gene_ids,
      reference_filter = annotation_dfs$current_gene_group_df$reference,
      filter_data = FALSE
      )
    
    print("merging")
    # merge
    msmt_wide <- interaction_models_wide %>% 
      full_join(singleplex_wide, by=gene_ids)
    
    print("filtering")
    
    # filter for current genes
    msmt_wide <- msmt_wide %>%
      inner_join(annotation_dfs$current_gene_group_df, by=c("reference"))
    
    # order to display nicely
    gene_id_cols <- c(gene_ids, "gene_group")
    other_cols <- setdiff(colnames(msmt_wide), gene_id_cols)
    msmt_wide <- msmt_wide[,c(gene_id_cols, other_cols)]
    
    return(msmt_wide)

  })
  
  ## Helper Functions
  # for gene inputs
  # 1. gene_group
  # 2. brush region
  # 3. gene symbol (in update_gene_inputs because simple)
  # plus function to update the reactive vals 
  load_gene_group_csv <- function(infile) {
    ext <- tools::file_ext(infile)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    # read the input
    new_gene_group_df <- read.csv(infile)
    if (ncol(new_gene_group_df) >= 1) {
      if (ncol(new_gene_group_df) >= 2) {
        new_gene_group_df <- new_gene_group_df[,1:2]
        colnames(new_gene_group_df) <- c("gene_id", "gene_group")
        
      } else {
        colnames(new_gene_group_df) <- c("gene_id")
        new_gene_group_df$gene_group <- "Group1"
      }
      
      new_gene_group_df$gene_id <- toupper(new_gene_group_df$gene_id)
    } else {
      new_gene_group_df <- data.frame(reference=character(0), gene_group=character(0))
    }
    
    return(new_gene_group_df)
    
  }
  
  get_genes_in_brush_region <- function(brushed_data) {
    
    brush_xcol <- brushed_data$mapping$x
    brush_ycol <- brushed_data$mapping$y

    
    brush_output_id <- brushed_data$outputId

    ds_panel <- gsub("region", "", strsplit(brush_output_id, "__")[[1]][1])
    dataset_type <- paste0(ds_panel, "_data")
    
    dataset_name <- strsplit(brush_output_id, "-")[[1]][3]
    # check if the dataset name is from a dropdown
    # if so get from input
    if (grepl("select", dataset_name)) {
      dataset_name <- input[[dataset_name]] 
    }

    if (!is.null(dataset_name)) {
      original_df <- get_dataset_from_type(dataset_type)[[dataset_name]]
      
      genes_in_region <- (original_df %>% 
                            dplyr::filter(.data[[brush_xcol]] < brushed_data$xmax &
                                            .data[[brush_xcol]] > brushed_data$xmin &
                                            .data[[brush_ycol]] < brushed_data$ymax & 
                                            .data[[brush_ycol]] > brushed_data$ymin))$reference
      
      
      new_gene_group_df <- data.frame(gene_id=genes_in_region, 
                                      gene_group="Regional_Selection")
    } else {
      new_gene_group_df <- data.frame(gene_id=character(0),
                                      gene_group=character(0))
    }
    
    return(new_gene_group_df)
    
  }
  
  # updates annotation_dfs/reactive vals
  update_gene_group_reactive_vals <- function(new_gene_group_df) {
    # updates the gene_group reactive vals
    # current_gene_group_df and gene_group_color_map
    # TODO: throw error if gene_symbol or reference is missing
    
    if (nrow(new_gene_group_df) > 0) {
    # now make gene color map
    unique_gene_groups <- unique(new_gene_group_df$gene_group)
    n_gene_groups <- length(unique_gene_groups)
    # add background in gray
    unique_gene_groups <- c("background", unique_gene_groups)
    gene_group_color_map <- accent_colors[1:(n_gene_groups+1)]
    names(gene_group_color_map) <- unique_gene_groups
    
    # update values
    annotation_dfs$current_gene_group_df <- new_gene_group_df
    annotation_dfs$gene_group_color_map <- gene_group_color_map
    } else {
      annotation_dfs$current_gene_group_df <- data.frame(reference=character(0),
                                                         gene_group=character(0))
      annotation_dfs$gene_group_color_map <- NULL
    }
    
    print("current annotation df")
    print(annotation_dfs$current_gene_group_df)
  }
  
  update_gene_input <- function(input_value, input_name, filter_bool=FALSE) {
    # checks gene_id_df with column gene_symbol against reference to determine if
    # gene name, uniprot, or ensembl
    # returns a list of the best id, intersection size, and 
    # df of updated ids and old id

    
    # input_name
    if (grepl("symbol", input_name)) {
      if (nchar(input_value) > 0) {
        
        # specify new gene group df
        gene_ids <- toupper(c(gsub(" ", "", input_value)))
        new_gene_group_df <- data.frame(gene_id=gene_ids,
                                        gene_group="Group1")
      } else {
        new_gene_group_df <- data.frame(gene_id=character(0),
                                        gene_group=character(0))
      }

    } else if (grepl("region", input_name)) {
      # expects gene_id, gene_group 
      new_gene_group_df <- get_genes_in_brush_region(input_value)
      
    } else if (grepl("group", input_name)) {
      file <- input_value
      # expects columns gene_id, gene_group
      new_gene_group_df <- load_gene_group_csv(file$datapath)
      
    }
    
      new_gene_group_df <- new_gene_group_df %>%
      subset(!is.na(gene_id) & length(gsub(" ", "", gene_id)) > 0)
    
    if (nrow(new_gene_group_df) > 0) {
      id_reference_columns <- c("reference", "gene_symbol", "uniprot", "ensembl")
      gene_ids <- new_gene_group_df$gene_id
      unique_intersection_size <- sapply(
        id_reference_columns, 
        function(col) length(intersect(as.character(gene_ref_map[[col]]),
                                       gene_ids)))

      best_idx <- unlist(which.max(unique_intersection_size))
      best_id <- id_reference_columns[best_idx]
      intersection_size <- unique_intersection_size[best_idx]

      
      # inner join will drop rows if missing!!
      # but interseciton size more meaningful and already reported
      # inner_join failed here for some unknown reason
      updated_gene_group_df <- merge(gene_ref_map, new_gene_group_df,
                                     by.x = best_id, by.y = "gene_id")
 
      new_gene_group_df <- updated_gene_group_df[,c("reference", "gene_group")]
      
      # remove nas in reference
      new_gene_group_df <- new_gene_group_df %>%
        subset(!is.na(reference))
      
      # if there are duplicated references, need to collapse here
      new_gene_group_df <- new_gene_group_df %>%
        group_by(reference) %>%
        summarize(gene_group=paste0(sort(unique(gene_group)), collapse=" & "))
      
      
    } else {
      new_gene_group_df <- data.frame(reference=character(0),
                                      gene_group=character(0))
      
    }
    

    update_gene_group_reactive_vals(new_gene_group_df)
    
    
    # update the apex_sp_data dfs
    reactive_df_name <- paste0(strsplit(input_name, "__")[[1]][1],
                               "_data")
    annotation_dfs[[reactive_df_name]] <- update_data_annotations(
      get_dataset_from_type(reactive_df_name),
      new_gene_annot_df=annotation_dfs$current_gene_group_df,
      filter_data=filter_bool)
    
  }


  ######################
  ## Main inputs
  ######################
  
  # Selection Inputs 
  # apex comparison
  apex_sp_comparison <- reactive({
    gsub(" vs ", "vs", input$apex_sp_comparison)
  })
  
  # update selectInput for apex time
  observeEvent(input$apex_time_dataset_select, {
    updateSelectInput(session, "apex_time_comparison",
                      choices = get_apex_time_comparison_choices(),
                      selected = get_apex_time_comparison_choices()[1]
    )
  })
  
  apex_time_dataset <- reactive({
    input$apex_time_dataset_select
  })
  
  apex_time_comparison <- reactive({
    input$apex_time_comparison
  })
  
  observeEvent(input$ip_filter, ignoreNULL = TRUE, {
    # need to update all the gene inputs for IP
    filter_status <- input$ip_filter == "Yes"
    suffixes <- c("symbol", "panelregion_plot_brush", "gene_group_input") 
    for (suffix in suffixes) {
      input_name <- paste0("ip__", suffix)
      input_val <- input[[input_name]]
      if (!is.null(input_val)) {
        update_gene_input(input_val, input_name, filter_bool=filter_status)
      }
    }
  })
  
  # Gene Inputs
  # observe gene symbols
  observeEvent(input$ip__symbol, ignoreNULL = TRUE, {
    filter_status <- input$ip_filter == "Yes"
    update_gene_input(input$ip__symbol, "ip__symbol", filter_bool=filter_status)
  })
  
  observeEvent(input$apex_sp__symbol, ignoreNULL = TRUE, {
    update_gene_input(input$apex_sp__symbol, "apex_sp__symbol")
    })
  
  observeEvent(input$apex_time__symbol, ignoreNULL = TRUE, {
    update_gene_input(input$apex_time__symbol, "apex_time__symbol")
  })
  
  # observe region brushes
  observeEvent(input$ip__panelregion_plot_brush, ignoreNULL = TRUE, {
    filter_status <- input$ip_filter == "Yes"
    update_gene_input(input$ip__panelregion_plot_brush, "ip__panelregion_plot_brush",
                      filter_bool = filter_status)
  })
  
  observeEvent(input$apex_sp__panelregion_plot_brush, ignoreNULL = TRUE, {
    update_gene_input(input$apex_sp__panelregion_plot_brush, "apex_sp__panelregion_plot_brush")
  })
  
  observeEvent(input$apex_time__panelregion_plot_brush, ignoreNULL = TRUE, {
    update_gene_input(input$apex_time__panelregion_plot_brush, "apex_time__panelregion_plot_brush")
  })
  
  # observe gene groups
  observeEvent(input$ip__gene_group_input, ignoreNULL=TRUE, {
    filter_status <- input$ip_filter == "Yes"
    update_gene_input(input$ip__gene_group_input, "ip__gene_group_input",
                      filter_bool = filter_status)
  })
  
  observeEvent(input$apex_sp__gene_group_input, ignoreNULL=TRUE, {
    update_gene_input(input$apex_sp__gene_group_input, "apex_sp__gene_group_input")
  })
  
  observeEvent(input$apex_time__gene_group_input, ignoreNULL=TRUE, {
    update_gene_input(input$apex_time__gene_group_input, "apex_time__gene_group_input")
  })
  
  
  ######################
  ## Outputs
  ######################
  
  ## Scatterplot outputs
  observe({
    panel_name <- get_current_panel()
    
      # name determines reactive or static scatter
      plot_type <- ifelse(grepl("region", panel_name), 
                          "scatterplot_reactive", "scatterplot_static")
      treat <- "all"
      
      if (grepl("ip", panel_name)) {
        ds_comp_grid <- expand.grid(
          dataset_name=unique(get_dataset_from_type("ip_data")$design$dataset_name),
          plotname_suffix=unique(get_dataset_from_type("ip_data")$comparisons$comparison)
        )
        ds_comp_grid$xcol <- paste0("log2FC_", ds_comp_grid$plotname_suffix)
        ds_comp_grid$ycol <- paste0("negLog10qval_", ds_comp_grid$plotname_suffix)
        ds_comp_grid$dataset_output_name <- ds_comp_grid$dataset_name
        
        ds_type <- "ip_data"
        
      } else if (grepl("apex_sp", panel_name)) {
        # define for single plex
        sp_comp_grid <- data.frame(dataset_name="all_bait_treat_single_plex",
                                   comparison=apex_sp_comparison())
        sp_comp_grid$xcol <- paste0("log2FC_", sp_comp_grid$comparison)
        sp_comp_grid$ycol <- paste0("negLog10qval_", sp_comp_grid$comparison)

        # define for linear models
        coef_names <- c("cellLine.CL097", "cellLine.Nigericin", "CL097.Nigericin")
        lm_comp_grid <- data.frame(
          comparison=coef_names,
          dataset_name=paste0(coef_names, "_linear_model"),
          xcol=paste0("b_", coef_names),
          ycol=paste0("negLog10qval_", coef_names)
        )
        
        ds_comp_grid <- rbind(sp_comp_grid, lm_comp_grid)   
        ds_comp_grid$plotname_suffix <- "any"
        ds_comp_grid$dataset_output_name <- ds_comp_grid$dataset_name
        
        ds_type <- "apex_sp_data"
        
      } else if (grepl("apex_time", panel_name)) {

        ds_comp_grid <- data.frame(
          dataset_name=apex_time_dataset(), 
          comparison=apex_time_comparison()
        )
        ds_comp_grid$xcol <- paste0("log2FC_", ds_comp_grid$comparison)
        ds_comp_grid$ycol <- paste0("negLog10qval_", ds_comp_grid$comparison)
        ds_comp_grid$dataset_output_name <- "apex_time_dataset_select"
        ds_comp_grid$plotname_suffix <- "any"
        
        ds_type <- "apex_time_data"
      }
      print(ds_comp_grid)
      lapply(1:nrow(ds_comp_grid),
             function(i) {
               print(ds_comp_grid$dataset_output_name[i])
               pltname <-
                 sprintf("%s-%s-%s-%s",
                         panel_name, plot_type, 
                         ds_comp_grid$dataset_output_name[i], 
                         ds_comp_grid$plotname_suffix[i])
               print(pltname)
               
               output[[pltname]] <-
                 renderPlot({
                   # if symbol or region plot gene symbols
                     plot_scatter(
                       annotation_dfs[[ds_type]],
                       ds_comp_grid$dataset_name[i],
                       get_dataset_from_type(ds_type)$design,
                       ds_comp_grid$xcol[i],
                       ds_comp_grid$ycol[i],
                       #annotate_genes = FALSE,
                       gene_group_color_map = annotation_dfs$gene_group_color_map
                       )
                   })
               })
  })
  
  ## Heatmaps outputs
  observe({
    panel_name <- get_current_panel()
    if (grepl("group|region", panel_name)) {
      
      plot_type <- "heatmap_static"
      treat <- "all"
      
      if (grepl("ip", panel_name)) {
        dataset_names <- 
          unique(get_dataset_from_type("ip_data")$design$dataset_name)
        ds_type <- "ip_data"
      } else if (grepl("apex_sp", panel_name)) {
        dataset_names <- c("all_bait_treat_single_plex")
        ds_type <- "apex_sp_data"
      } else if (grepl("apex_time", panel_name)) {
        dataset_names <- c("nlrp3_nigericin_timecourse", "nlrp3_cl097_timecourse",
                           "p4c_nigericin_timecourse", "p4c_cl097_timecourse")
        ds_type <- "apex_time_data"
      }
      
      lapply(dataset_names,
             function(ds_name) {
               pltname <-
                 sprintf("%s_%s_%s_%s",
                         panel_name, plot_type, ds_name, treat)
               output[[pltname]] <-
                 renderPlot({
                   plot_heatmap(annotation_dfs[[ds_type]],
                            ds_name,
                            get_dataset_from_type(ds_type)$design,
                            gene_group_color_map = annotation_dfs$gene_group_color_map
               )
               
             })
             })
    }
    })
  
  ## Boxplot outputs
  observe({
    panel_name <- get_current_panel()
    if (grepl("symbol", panel_name) & grepl("ip|apex_sp", panel_name)) {
      
      plot_type <- "boxplot_static"
      treat <- "all"
      
      if (grepl("ip", panel_name)) {
        dataset_names <- 
          unique(get_dataset_from_type("ip_data")$design$dataset_name)
        ds_type <- "ip_data"
      } else if (grepl("apex_sp", panel_name)) {
        dataset_names <- c("all_bait_treat_single_plex")
        ds_type <- "apex_sp_data"
      } 
      
      lapply(dataset_names,
             function(ds_name) {
               pltname <-
                  sprintf("%s-%s-%s-%s",
                          panel_name, plot_type, ds_name, treat) 

               output[[pltname]] <-
                 renderPlot({
                   plot_boxplot_intensity(
                     annotation_dfs[[ds_type]],
                     ds_name,
                     get_dataset_from_type(ds_type)$design,
                     gene_group_color_map = annotation_dfs$gene_group_color_map
                   )
                   
                 })
             })
    }
  })
  

  ## Violinplot outputs
  observe({
    panel_name <- get_current_panel()
    plot_type <- "violinplot_static"
    
    if (grepl("region|group", panel_name) & grepl("ip|sp", panel_name)) {
    plotname_suffix <- "all"
    
    if (grepl("ip", panel_name)) {
      ds_comp_grid <- expand.grid(
        dataset_name=unique(get_dataset_from_type("ip_data")$design$dataset_name),
        msmt_col_regex=c("log2FC_"),
        ylabel=c("Log2 Fold Change")
      )
      ds_comp_grid$dataset_output_name <- ds_comp_grid$dataset_name
      
      ds_type <- "ip_data"
      
    } else if (grepl("apex_sp", panel_name)) {
      # define for single plex
      sp_comp_grid <- data.frame(dataset_name="all_bait_treat_single_plex",
                                 msmt_col_regex=paste0("log2FC_", apex_sp_comparison()),
                                 ylabel=("Log2 Fold Change"))
      
      # define for linear models
      coef_names <- c("cellLine.CL097", "cellLine.Nigericin", "CL097.Nigericin")
      lm_comp_grid <- data.frame(
        dataset_name=paste0(coef_names, "_linear_model"),
        msmt_col_regex=paste0(paste0("b_", coef_names), collapse="|"),
        ylabel=c("Log2 Coefficient")
      )
      
      ds_comp_grid <- rbind(sp_comp_grid, lm_comp_grid)   
      ds_comp_grid$dataset_output_name <- ds_comp_grid$dataset_name
      
      ds_type <- "apex_sp_data"
      
    }
    
    lapply(1:nrow(ds_comp_grid),
           function(i) {
             pltname <-
               sprintf("%s-%s-%s-%s",
                       panel_name, plot_type, 
                       ds_comp_grid$dataset_output_name[i], 
                       plotname_suffix)
             
             # get the measurement columns
             msmt_cols <- get_measurement_columns(
               annotation_dfs[[ds_type]][[ds_comp_grid$dataset_name[i]]],
               ds_comp_grid$msmt_col_regex)

             output[[pltname]] <-
               renderPlot({
                 # if symbol or region plot gene symbols
                 plot_violin_msmt(
                   annotation_dfs[[ds_type]],
                   ds_comp_grid$dataset_name[i],
                   msmt_cols,
                   gene_group_color_map = annotation_dfs$gene_group_color_map
                 ) + 
                   geom_hline(yintercept=0, color='black', linetype='dashed') +
                   ylab(ds_comp_grid$ylabel[i])

               })
           })
    }
  })
  

  ## Timecourse lm plots
  observe({
    panel_name <- get_current_panel()
    if (grepl("time", panel_name)) {
      dataset_names <- c("nigericin_timecourse_linear_model", "cl097_timecourse_linear_model")
      plot_type <- ifelse(grepl("region", panel_name), "scatterplot_reactive", "scatterplot_static")

      lapply(dataset_names,
             function(ds_name) {
               pltname <- sprintf("%s-%s-%s-any", panel_name, plot_type, ds_name)
               output[[pltname]] <-
                 renderPlot({
                   lm_df <- annotation_dfs$apex_time_data[[ds_name]]
                   plot_timecourse_linear_model(
                     lm_df,
                     gene_group_color_map = annotation_dfs$gene_group_color_map
                   )
                 })
             })
      
    }
  })

  ## Data table output
  observe({
    panel_name <- get_current_panel()
    tbl_name <- paste0(panel_name, "_table")
    
    output[[tbl_name]] <-
      DT::renderDataTable({
        
        if (grepl("ip", panel_name)) {
          df <- ip_data_table_output()
        } else if (grepl("apex_sp", panel_name)) {
          df <- apex_sp_data_table_output()
        } else if (grepl("apex_time", panel_name)) {
          df <- apex_time_data_table_output()
        }
        
        # set the reactive value
        annotation_dfs$data_table_output <- df
        
        DT::datatable(df,
                      extensions = 'FixedColumns',
                      options =
                        list(scrollX = T,
                             fixedColumns = list(leftColumns = 3))
      )
      })
  })
  
  ######################
  ## Data download
  ######################
  
  output$downloadData <- downloadHandler(
      filename = "gene_group_statistics.csv",
      content = function(file) {
        write.csv(annotation_dfs$data_table_output,
                  file,
                  row.names = FALSE
        )
      }
    )
        
}

# Run the application
shinyApp(ui = ui, server = server, options = list(width = '100%'))