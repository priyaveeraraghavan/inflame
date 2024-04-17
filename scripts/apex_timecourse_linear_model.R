library(dplyr)

working_directory <- "/Users/priyaveeraraghavan/Documents/GitHub/inflame/paper_analysis/"

# define the bait-combined normalized files
combined_normalized_files <- data.frame(
  filename=file.path(working_directory, "data",
                     c("all_stats_results_P4CNL3_combined_Nig_1peptide_medianpolish.csv",
                       "all_stats_results_NL3P4C_combined_CL097_1peptide_medianPolish.csv")),
  treatment=c("Nigericin", "CL097"))

# load the bait-combined normalized dfs
combined_normalized_files$dfs <- lapply(1:nrow(combined_normalized_files),
                                        function(x) read.csv(combined_normalized_files$filename[x], 
                                                             stringsAsFactors = FALSE) %>%
                                          rowwise() %>%
                                          mutate(treatment=combined_normalized_files$treatment[x]))



# create the design matrix for different conditions (CTR vs NLRP)
get_design_matrix_condition <- function(input_data) {
  # Gets the design information by parsing colnames of input_data
  # Arguments: 
  #  input_data: data.frame with MSStats dataProcess output intensities, 
  #              columns must be the rownames of design_mat
  
  design_mat <- data.frame(
    sample_name=colnames(input_data)[grepl("^CTR|^NLR", 
                                           colnames(input_data))])
  
  design_mat$condition <- gsub("[_0-9]+", "", design_mat$sample_name)
  design_mat$replicate <- sapply(design_mat$sample_name, 
                                 function(x) strsplit(x, "_")[[1]][2])
  design_mat <- data.frame(subset(design_mat, !is.na(replicate)))
  design_mat$timepoint <- sapply(design_mat$sample_name, 
                                 function(x) as.numeric(strsplit(gsub("CTR|NLRP", "", x), "_")[[1]][1]))
  rownames(design_mat) <- design_mat$sample_name
  design_mat$sample_name <- NULL
  
  return(design_mat)
}


# get the intensity mat
get_intensity_matrix <- function(input_data, design_mat) {
  # Gets the intensity matrix from input_data and desired columns (from design_mat)
  # Arguments: 
  #  input_data: data.frame with MSStats dataProcess output intensities, 
  #              columns must be the rownames of design_mat
  #  design_mat: data.frame design matrix with sample names as rownames. Must
  #              be the same rownames as the columns of input_data
  #
  # Returns: intensity data.frame of size n_genes x n_samples, where colnames
  #          are the sample names and rownames are the reference gene names
  
  intensity_data <- data.frame(input_data[,rownames(design_mat)])
  rownames(intensity_data) <- input_data$Reference
  
  return(intensity_data)
}


combined_normalized_files$design_mat <- lapply(combined_normalized_files$dfs,
                                               function(x) get_design_matrix_condition(x) %>%
                                                 subset(timepoint <= 30))

combined_normalized_files$intensity_mat <- lapply(1:nrow(combined_normalized_files),
                                                  function(x) get_intensity_matrix(combined_normalized_files$dfs[x][[1]],
                                                                                   combined_normalized_files$design_mat[x][[1]]))


get_valid_genes <- function(intensity_data,
                            design,
                            num_reps_nonzero=3) {
  design$sample_name <- rownames(design)
  intensity_data$reference_name <- rownames(intensity_data)
  df_long <- intensity_data %>% 
    pivot_longer(contains("CTR") | contains("NLRP") | contains("t"), 
                 names_to = "sample_name", values_to = "intensity")
  
  df_long <- df_long %>% inner_join(design, by="sample_name")
  
  valid_condition <- df_long %>%
    group_by(reference_name, condition, timepoint) %>%
    summarize(is_valid_condition=sum(intensity > 0) >= num_reps_nonzero)
  
  valid_gene <- valid_condition %>%
    group_by(reference_name) %>%
    summarize(is_valid=sum(!is_valid_condition) == 0,
              is_valid_nlrp3=sum(!is_valid_condition[grepl("NLRP", condition)]) == 0,
              is_valid_control=sum(!is_valid_condition[grepl("CTR", condition)]) == 0)
  
  return(valid_gene)
}


combined_normalized_files$valid_genes <- lapply(1:nrow(combined_normalized_files),
                                                function(x) get_valid_genes(combined_normalized_files$intensity_mat[x][[1]],
                                                                            combined_normalized_files$design_mat[x][[1]]))



# define fitting function
fit_single_protein <- function(reference_name, 
                               intensity_data, 
                               design,
                               formula_string="abundance ~ condition + timepoint + condition*timepoint",
                               return_lm=FALSE,
                               to_plot=FALSE) {
  
  # TODO: check if reference is present
  # TODO: check that the intensity data sample names are the same as in the design
  intensity_data_gene_long <- intensity_data[reference_name,]
  intensity_data_gene <- data.frame(t(intensity_data_gene_long))
  colnames(intensity_data_gene) <- "abundance"
  
  model_df <- merge(intensity_data_gene, design, by='row.names', all.y=TRUE)
  model_df$timepoint <- model_df$timepoint #/10
  full_fit <- lm(as.formula(formula_string), 
                 data=model_df)
  
  model_df$prediction <- full_fit$fitted.values
  
  if (to_plot==TRUE) {
    
    plot(ggplot(model_df, aes(timepoint, abundance, color=condition)) + 
           geom_point())
  }
  
  # process the fit summary
  coefs <- data.frame(summary(full_fit)$coefficients)
  colnames(coefs) <- c("estimate", "stderr", "tvalue", "pvalue")
  coefs$beta <- gsub("\\(|\\)", "", rownames(coefs))
  coefs_wide <- data.frame(coefs %>% 
                             pivot_wider(names_from=beta, 
                                         values_from=c(estimate, stderr, tvalue, pvalue)))
  coefs_wide$reference <- reference_name
  if (return_lm) {
    coefs_wide$lm_full_fit <- list(full_fit)
  }
  
  return(list(model_df=model_df, coefs=coefs_wide, model_fit=full_fit))
  
}

# run timecourse lms
combined_normalized_files <- combined_normalized_files %>%
  rowwise() %>%
  mutate(interaction_model=list(data.frame(
    do.call(rbind, 
            lapply(subset(valid_genes, is_valid)$reference_name,
                   function(ref) fit_single_protein(ref, 
                                                    intensity_mat,
                                                    design_mat,
                                                    formula_string="abundance ~ timepoint + condition + timepoint*condition")$coefs
            )
    )
  )
  )
  )

# get adjusted pvalues
calculate_adjusted_pvalues <- function(model_results) {
  
  pval_colnames <- colnames(model_results)[grepl("pvalue", colnames(model_results))]
  betas <- gsub("pvalue_", "", pval_colnames)
  for (beta in betas) {
    pval_colname <- sprintf("pvalue_%s", beta)
    padj_colname <- sprintf("qvalue_%s", beta)
    model_results[padj_colname] <- p.adjust(
      unlist(model_results[pval_colname]), method="BH"
    )
  }
  
  return(model_results)
  
}

combined_normalized_files$interaction_model <- lapply(combined_normalized_files$interaction_model,
                                                      calculate_adjusted_pvalues) 


# merge output data and write to file
for (i in c(1,2)) {
  
  # add "intensity" to colnames
  int_mat <- combined_normalized_files$intensity_mat[[i]]
  colnames(int_mat) <- gsub("_", "_Intensity_", colnames(int_mat))
  # make reference a column for merging
  int_mat$reference <- rownames(int_mat)
  
  # get the gene names
  int_mat$gene_symbol <- combined_normalized_files$dfs[[i]]$Gene.Symbol
  
  # get min and max times used by model. Useful for plotting later
  designs <-  combined_normalized_files$design_mat[[i]]
  min_time <- min(designs$timepoint)
  max_time <- max(designs$timepoint)
  
  outdf <- int_mat %>%
    inner_join(combined_normalized_files$interaction_model[[i]], by='reference')
  
  outdf$min_time <- min_time
  outdf$max_time <- max_time
  treat_name <- combined_normalized_files$treatment[[i]]
  
  write.csv(outdf, 
            file=file.path(working_directory, "data", 
                           sprintf("%s_timecourse_linear_model.csv", treat_name)))
  
  
}

                                                  