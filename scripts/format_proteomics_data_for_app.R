# script to build aggregated organelle IP data
# broadly
# loads the filtered and unfiltered IP datasets
# annotates the unfiltered dataset with whether a gene is found in the filtered data
# annotates the unfiltered dataset with the organelle
# joins all the unfiltered datasets with the appropriate annotations
# saves for future loading

library(tidyr)
library(dplyr)
library(ggrastr)
library(ggplot2)

project_dir <- "/Users/priyaveeraraghavan/Documents/GitHub/inflame/Proteomics_Visualization"

data_dir <- file.path(project_dir, "data")

###########
## Helper/universal
########

clean_gene_names <- function(df) {

  colnames(df) <- gsub("Gene.Symbol", "gene_symbol", colnames(df))
  colnames(df) <- gsub("Reference", "reference", colnames(df))
  
  # remove SP from reference
  df$reference <- gsub("sp\\|", "", df$reference)
  
  # move gene_symbol to gene_symbol_species
  df$gene_symbol_species <- df$gene_symbol
  
  # create an uppercase universal gene symbol
  df$gene_symbol <- toupper(df$gene_symbol)
  
  return(df)
  
}

############ 
#### IP Data
############

ip_metadata <- read.csv(file.path(data_dir, "IP_metadata.csv")) 


load_annotate_organelle <- function(organelle, metadata=ip_metadata) {
  
  # Subset metadata to include only rows with the specified organelle
  organelle_df <- subset(metadata, dataset_name == organelle)
  
  # Check if the organelle is "whole_cell"
  if (organelle == "whole_cell") {
    # Subset metadata for "whole_cell"
    meta_df <- subset(metadata, dataset_name == organelle)
    # Construct file path using metadata information
    fname <- file.path(data_dir, meta_df$folder, meta_df$filename)
    
    # Read CSV file into prot_df
    prot_df <- read.csv(fname)
    # Add columns to prot_df
    prot_df$dataset_name <- organelle
    prot_df$filtered <- TRUE
  } else {
    # Process for other organelles
    
    # Unfiltered data processing
    unfiltered_meta_df <- subset(metadata, dataset_name == organelle & grepl("unfiltered", folder))
    unfiltered_fname <- file.path(data_dir, unfiltered_meta_df$folder, unfiltered_meta_df$filename)
    print(unfiltered_fname)
    unfiltered_IP_df <- read.csv(unfiltered_fname)
    unfiltered_IP_df$dataset_name <- organelle
    
    # Filtered data processing
    filtered_meta_df <- subset(metadata, dataset_name == organelle & !grepl("unfiltered", folder))
    filtered_fname <- file.path(data_dir, filtered_meta_df$folder, filtered_meta_df$filename)
    filtered_IP_df <- read.csv(filtered_fname)
    print(filtered_fname)
    
    # Add if filtered or not to unfiltered df
    unfiltered_IP_df$filtered <- unfiltered_IP_df$Reference %in% filtered_IP_df$Reference
    prot_df <- unfiltered_IP_df
  }
  
  # clean prot df columns
  prot_df <- clean_gene_names(prot_df)

  # Modify column names
  colnames(prot_df) <- gsub("_vs_", "vs", colnames(prot_df))
  
  # Compute negative log10 values for qval columns
  for (qval in colnames(prot_df)[grepl("qval", colnames(prot_df))]) {
    prot_df[,sprintf("negLog10%s", qval)] <- -log10(prot_df[,qval])
  }
  
  # Print uncharacterized fragments
  print(prot_df %>% subset(nchar(gene_symbol) == 0))
  # Filter out "#" in the reference column and empty gene symbols
  prot_df <- prot_df %>% subset(!grepl("#", reference) & nchar(gene_symbol) > 0) 
  
  # Print information about data validity
  prev_shape <- nrow(prot_df)
  new_shape <- nrow(prot_df)
  print(sprintf("%s out of %s annotations valid", new_shape, prev_shape))
  
  # Return the processed dataframe
  return(prot_df)
}


ip_datasets <- unique(ip_metadata$dataset_name)
processed_ip_data <- lapply(ip_datasets, load_annotate_organelle)
names(processed_ip_data) <- ip_datasets

# process the metadata to get design info 
# columns: compartment / condition / replicate 
get_design <- function(ds_name, processed_data) {
  dataset_df <- processed_data[[ds_name]] 
  
  # intensity cols contain the design info
  intensity_cols <- colnames(dataset_df)[grepl("Intensity", colnames(dataset_df))] 
  
  # conditions will be Control, Nigericin, and CL097
  conditions <- sapply(intensity_cols, function(x) strsplit(x, "_Intensity_")[[1]][1])
  
  # replicates will be 1, 2, 3, etc and Mean
  replicates <- sapply(intensity_cols, function(x) strsplit(x, "_Intensity_")[[1]][2])
  

  res <- data.frame(dataset_name=ds_name,
                    condition=conditions,
                    replicate=replicates,
                    column_name=intensity_cols)
  
  res <- res %>%
    rowwise() %>%
    mutate(row_names=sprintf("%s_%s_%s", ds_name, condition, replicate))
  
  # cast as df
  res <- data.frame(res)
  
  rownames(res) <- res$row_names
  res$row_names <- NULL
  
  return(res)
}

design_info <- do.call(rbind, lapply(ip_datasets, 
                                     function(ds) get_design(ds, processed_ip_data)))
processed_ip_data[['design']] <- design_info

# get comparisons info
get_comparisons <- function(ds_name, processed_data) {
  dataset_df <- processed_data[[ds_name]] 
  
  fc_cols <- colnames(dataset_df)[grepl("log2FC", colnames(dataset_df))]
  
  # comparisons
  comparisons=gsub("log2FC_", "", fc_cols)
  
  # just for clarity call out numerator and denominator
  numerators <- sapply(comparisons, function(x) strsplit(x, "vs")[[1]][1])
  denominators <- sapply(comparisons, function(x) strsplit(x, "vs")[[1]][2])
  
  res <- data.frame(dataset_name=ds_name,
                    comparison=comparisons, 
                    numerator=numerators,
                    denominator=denominators)
  
  rownames(res) <- NULL
  
  return(res)
  
  
}

comparison_info <- do.call(rbind, lapply(ip_datasets, 
                                         function(ds) get_comparisons(ds, processed_ip_data)))
processed_ip_data[['comparisons']] <- comparison_info


## Make pre-rendered scatterplots
rasterize_scatter <- function(compartment, comparison, filter_data) {
  ###
  xcol <- paste0(c("log2FC", comparison), collapse="_")
  ycol <- paste0(c("negLog10qval", comparison), collapse="_")
  
  data <- processed_ip_data[[compartment]] 
  print(compartment)
  print(data)
  if (filter_data) {
    data <- subset(data, filtered)
  }
  
  plt <- rasterize(ggplot(data, aes_string(xcol, ycol)) +
    geom_point(show.legend=FALSE), dpi=100) 
  
  return(plt)
  
}

premake_rasterized_scatters <- function(compartment) {
  comparisons <- unique(processed_ip_data[['comparisons']]$comparison)
  filtered_scatter <- lapply(comparisons, function(compare) rasterize_scatter(compartment, compare, TRUE))
  names(filtered_scatter) <- comparisons
  unfiltered_scatter <- lapply(comparisons, function(compare) rasterize_scatter(compartment, compare, FALSE))
  names(unfiltered_scatter) <- comparisons
  return(list(filtered=filtered_scatter, unfiltered=unfiltered_scatter))
}

#compartments <- unique(processed_ip_data$design$compartment)
#processed_ip_data[['rasterized_scatter']] <- lapply(compartments, premake_rasterized_scatters)
#names(processed_ip_data[['rasterized_scatter']]) <- compartments
  
################### 
#### APEX Data ####
###################
apex_metadata <- read.csv(file.path(data_dir, "APEX_metadata.csv"))

#### Timecourse ###
# Metadata
apex_timecourse_metadata <- subset(apex_metadata, dataset_type == "timecourse" & bait != "all")

# Load and format the data
load_annotate_apex <- function(ds_name, metadata) {
  
  # Subset metadata to include only rows with the specified bait and treatment
  meta_subs <- subset(metadata, dataset_name == ds_name)
  
  prot_df <- read.csv(file.path(data_dir, meta_subs$folder, meta_subs$filename))

  
  # Rename columns in prot_df
  prot_df <- clean_gene_names(prot_df)

  # add info about dataset
  prot_df$dataset_name <- ds_name
  prot_df$filtered <- TRUE # by default is true here, for compatibility
  
  # Compute negative log10 values for qval columns
  for (qval in colnames(prot_df)[grepl("qval", colnames(prot_df))]) {
    prot_df[,sprintf("negLog10%s", qval)] <- -log10(prot_df[,qval])
  }
  
  # change - to vs
  colnames(prot_df) <- gsub("_vs_", "vs", colnames(prot_df))
  
  # Print uncharacterized fragments
  print(prot_df %>% subset(nchar(gene_symbol) == 0))
  # Filter out "#" in the reference column and empty gene symbols
  prev_shape <- nrow(prot_df)
  prot_df <- prot_df %>% subset(!grepl("#", reference) & nchar(gene_symbol) > 0) 
  
  # Print information about data validity
  new_shape <- nrow(prot_df)
  print(sprintf("%s out of %s annotations valid", new_shape, prev_shape))
  
  # Return the processed dataframe
  return(prot_df)
}

time_datasets <- apex_timecourse_metadata$dataset_name
processed_apex_timecourse_data <- lapply(time_datasets, 
                                         function(ds) load_annotate_apex(ds, apex_timecourse_metadata))
names(processed_apex_timecourse_data) <- time_datasets

# Design info
apex_time_design_info <- do.call(rbind, lapply(time_datasets, 
                                     function(ds) get_design(ds, processed_apex_timecourse_data)))

processed_apex_timecourse_data[['design']] <- apex_time_design_info

# Comparisons info
apex_time_comparison_info <- do.call(rbind, lapply(time_datasets, 
                                         function(ds) get_comparisons(ds, processed_apex_timecourse_data)))
processed_apex_timecourse_data[['comparisons']] <- apex_time_comparison_info

# grab the linear model info
apex_timecourse_metadata_joint <- subset(apex_metadata, dataset_type == "timecourse" & bait == "all")

for (ds_name in  unique(apex_timecourse_metadata_joint$dataset_name)) {
  ds_info <- subset(apex_timecourse_metadata_joint, dataset_name == ds_name)
  
  infile_name <- file.path(data_dir, ds_info$folder, ds_info$filename)
  
  prot_df <- read.csv(infile_name)
  
  # clean prot df columns
  prot_df <- clean_gene_names(prot_df)
  
  prot_df$filtered <- TRUE
  
  processed_apex_timecourse_data[[ds_name]] <- prot_df
}


#### Single-Plex ###
processed_apex_singleplex_data <- list()
processed_apex_singleplex_data[["all_bait_treat_single_plex"]] <- load_annotate_apex("all_bait_treat_single_plex", apex_metadata)

processed_apex_singleplex_data[["design"]] <- get_design("all_bait_treat_single_plex", processed_apex_singleplex_data)
processed_apex_singleplex_data[["comparisons"]] <- get_comparisons("all_bait_treat_single_plex", processed_apex_singleplex_data)

apex_sp_metadata_lm <- subset(apex_metadata, dataset_type == "single_plex" & grepl("linear_model", dataset_name))

for (ds_name in apex_sp_metadata_lm$dataset_name) {
  ds_info <- subset(apex_sp_metadata_lm, dataset_name == ds_name)
  
  infile_name <- file.path(data_dir, ds_info$folder, ds_info$filename)
  
  prot_df <- read.csv(infile_name) 
  
  # change reference column to have uniprot and descriptive
  prot_df <- prot_df %>%
    rowwise() %>%
    mutate(Reference=sprintf("%s|%s", Reference, GeneNameExtra))
  
  # fix column names
  prot_df <- clean_gene_names(prot_df)
  
  # add neglog10 qvals
  # Compute negative log10 values for qval columns
  colnames(prot_df) <- gsub("qvalue", "qval", colnames(prot_df))
  for (qval in colnames(prot_df)[grepl("qval", colnames(prot_df))]) {
    prot_df[,sprintf("negLog10%s", qval)] <- -log10(prot_df[,qval])
  }
  
  # set filtered to true
  prot_df$filtered <- TRUE
  processed_apex_singleplex_data[[ds_name]] <- prot_df
}

## get all gene refs
# get gene reference map
# includes all genes from all studies
# probably in the future want to annotate with more info
get_gene_reference_map <- function(ds_names, processed_data) {
  all_genes <- do.call(rbind, lapply(ds_names, 
                                     function(ds_name) {
                                       print(ds_name)
                                       df <- processed_data[[ds_name]][,c("gene_symbol", "reference", "gene_symbol_species")]
                                       print(head(df))
                                       print(df$reference[grepl("sp", df$reference)])
                                       df
                                       }))
  
  return(data.frame(unique(all_genes)))
}

gene_ref_map <- unique(data.frame(
  do.call(rbind, 
          lapply(list(processed_ip_data, processed_apex_singleplex_data, processed_apex_timecourse_data),
                 function(data) {
                   valid_ds <- names(data)[!grepl("design|comparisons", names(data))]
                   return(get_gene_reference_map(valid_ds, data))
                 }
          )
  )        
))

# add uniprot
# add ensembl
gene_ref_map$uniprot <- sapply(gene_ref_map$reference, function(x) strsplit(gsub("sp\\|", "", gsub("tr\\|", "", x)), "\\|")[[1]][1])
gene_ref_map$species <- sapply(gene_ref_map$reference, function(x) strsplit(x, "_")[[1]][2])
# update reference column
gene_ref_map <- gene_ref_map %>% dplyr::rename(old_reference=reference)
gene_ref_map <- gene_ref_map %>%
  rowwise() %>%
  mutate(reference=sprintf("%s|%s_%s", uniprot, gene_symbol, species))

mouse_mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
human_mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

ensembl_ids_mouse <- biomaRt::getBM(c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id"),
                                 filters="uniprot_gn_id", 
                                 values=subset(gene_ref_map, grepl("MOUSE", reference))$uniprot,
                                 mart=mouse_mart)

ensembl_ids_human <- biomaRt::getBM(c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id"),
                                    filters="uniprot_gn_id", 
                                    values=subset(gene_ref_map, grepl("HUMAN", reference))$uniprot,
                                    mart=human_mart)

ensembl_ids_combined <- rbind(ensembl_ids_mouse, ensembl_ids_human)
ensembl_ids_combined <- ensembl_ids_combined %>% dplyr::rename(ensembl=ensembl_gene_id, gene_symbol_species=external_gene_name)

gene_ref_map_updated <- gene_ref_map[,c("reference", "old_reference", "uniprot", "species")] %>% left_join(ensembl_ids_combined, by = c("uniprot" = "uniprot_gn_id"))
gene_ref_map_updated$gene_symbol <- toupper(gene_ref_map_updated$gene_symbol_species)

# remove nas
gene_ref_map_updated <- gene_ref_map_updated %>%
  subset(!is.na(reference))

# go back and fix the reference names in datasets
fix_ds_refs <- function(ds_obj_lst) {
  ds_names <- names(ds_obj_lst)[!grepl("comparisons|design", names(ds_obj_lst))]
  unchanged_names <- setdiff(names(ds_obj_lst), ds_names)
  new_obj_lst <- lapply(
    ds_names, 
    function(ds_name) {
      print(ds_name)
      df_original <- ds_obj_lst[[ds_name]]
      df_original$gene_symbol <- NULL
      print(dim(df_original))
      df_new <- gene_ref_map_updated[,c("reference", "old_reference")] %>% 
        # join based on old reference
        right_join(df_original, 
                  by=c("old_reference" = "reference")) %>%
        # remove old ref
        dplyr::select(-old_reference) %>%
        # drop any duplicates
        unique()
      print(dim(df_new))
      
      return(df_new)
                        })
  names(new_obj_lst) <- ds_names
  
  for (name in unchanged_names) {
    new_obj_lst[[name]] <- ds_obj_lst[[name]]
  }
  return(new_obj_lst)
}

fixed_processed_ip_data <- fix_ds_refs(processed_ip_data)
fixed_processed_apex_timecourse_data <- fix_ds_refs(processed_apex_timecourse_data)
fixed_processed_apex_singleplex_data <- fix_ds_refs(processed_apex_singleplex_data)


# now clean the gene_ref_map
gene_ref_map <- gene_ref_map_updated %>%
  dplyr::select(-old_reference) %>%
  unique()

save(fixed_processed_ip_data,
     fixed_processed_apex_timecourse_data, 
     fixed_processed_apex_singleplex_data,
     gene_ref_map,
     file=file.path(data_dir, "processed_proteomics_data.RData"))


