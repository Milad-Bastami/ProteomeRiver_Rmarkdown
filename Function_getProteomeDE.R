#---------- Proteome DE Analysis
#----- Packages
# suppressPackageStartupMessages({
#   p_load(magrittr)
#   p_load(limma)
#   p_load(qvalue)
#   p_load(ruv)
#   p_load(mixOmics)
#   p_load(ggrepel)
#   p_load(ggpubr)
#   p_load(statmod)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir = "results/proteomics/de_proteins"
# #- Minimum log2-fold-change below which changes not considered. Limma treat function.
# treat_lfc_cutoff = NA
# #- Normalization method, choices are none, scale, quantile or cyclicloess, default scale
# normalization = "cyclicloess"
# #- Mark zero values as NA (mark_as_na), or add a very small pseudo count (add_pseudo_count)
# handling_zeros = "mark_as_na"
# #- Implement imputation, logical
# imputation = FALSE
# #- A double from 0 to 1 indicating the percentage of quantified values in any treatment group
# impute_min_percent = NA
# #- An integer indicating the number of replicatses, in which the percent filter must be satisfied
# #- for the protein / phosphosite to be included
# impute_min_num_of_groups = NA
# #- A double from 0 to 1 indicating the minimum percentage of quantified values in any treatment
# #- group to be included in site and condition specific impuatation
# impute_specific_percent = NA
# #- logical, whether comparisons involving imputed values should be removed
# remove_imputed = FALSE
# #- A double from 0 to 1 indicating the minimum percentage of quantified values in any treatment
# #- group to be included in statistical pairwise comparison of groups
# remove_imputed_perc = 1
# #- logical, should an intensity-trend be allowed for the prior variance? Default is that the
# #- prior variance is constant
# eBayes_trend = TRUE
# #- logical, should the estimation of df.prior and var.prior be robustified against outlier
# #- sample variances?
# eBayes_robust = TRUE
# #- Remove protein if it exceeds this maximum number of samples with missing values per
# #- experimental group
# max_num_samples_miss_per_group = 0
# #- Abundance threshold above which the protein in the sample is accepted for analysis
# abundance_threshold = 0
# #- Regular expression pattern to identify columns with abundance values belonging to the
# #- experiment, default is "\\d+"
# group_pattern = "\\d+"
# #- q-value threshold below which a protein has statistically significant differetial expression
# q_val_thresh = 0.05
# #- q-value threshold below which a protein has statistically significant differetial expression,
# #- used for control genes
# control_genes_q_val_thresh = 0.05
# #- The number of unwanted factors to use
# ruv_k = 4
# #- The number of negative control proteins to use
# num_neg_ctrl = 100
# #- A string representing the ruv3 method to use
# ruv_method = "ruv3"
# #- Input file with the protein abundance values
# counts_table_file = "results/proteomics/clean_proteins/counts_table_cleaned.tab"
# #- Input file with a table listing all comparisons to be made in string, one comparison per line
# #- (e.g. groupB.vs.group_A = groupB - groupA)
# contrasts_file = "data/proteomics_contrast.tab"
# #- A string representing the formula for input into the model.frame function. (e.g. ~ 0 + group)
# formula_string = "~ 0 + group "
# #- Input file with the design matrix
# design_matrix_file = "data/proteomics_design.tab"
# #- A string describing the sample ID. This must be a column that exists in the design matrix
# sample_id = "Sample_ID"
# #- A string describing the replicate group ID. This must be a column that exists in the design matrix
# replicate_group_id = NA
# #- A string describing the experimental group ID. This must be a column that exists in the design matrix
# group_id = "group"
# #- A string describing the technical replicates that needs to be averaged together.
# #- This must be a column that exists in the design matrix
# average_replicates_id = NA
# #- A string describing the row id
# row_id = "uniprot_acc"
# #- A string to indicate the type of analysis and is used in the file name of the output results table
# file_prefix = "de_proteins"
# #- A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg]
# plots_format = list("pdf","png")

#----- Function
getProteomeDE <- function(
    output_dir,
    treat_lfc_cutoff,
    normalization,
    handling_zeros,
    imputation,
    impute_min_percent,
    impute_min_num_of_groups,
    impute_specific_percent,
    remove_imputed,
    remove_imputed_perc,
    eBayes_trend,
    eBayes_robust,
    max_num_samples_miss_per_group,
    abundance_threshold,
    group_pattern,
    q_val_thresh,
    control_genes_q_val_thresh,
    ruv_k,
    num_neg_ctrl,
    ruv_method,
    counts_table_file,
    contrasts_file,
    formula_string,
    design_matrix_file,
    sample_id,
    replicate_group_id,
    group_id,
    average_replicates_id,
    row_id,
    file_prefix,
    plots_format) {
  #----- Create Results Directories
  createDirectoryIfNotExists <- function(file_path, mode = "0777") {
    #- Create directory recursively if it doesn't exist
    if (!file.exists(file_path)) {
      dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)
    }
  }
  createDirectoryIfNotExists(file_path = output_dir)

  #----- Add average_replicates_id
  #- average_replicates_id <- "group"
  if (is.na(replicate_group_id) & (!is.na(group_id) | group_id == "")) {
    replicate_group_id <- group_id
  }

  #----- Contrasts File
  contrasts_tbl <- vroom::vroom(contrasts_file, delim = "\t")

  #----- Matrix File
  design_mat_cln <- vroom::vroom(design_matrix_file) %>%
    as.data.frame() %>%
    dplyr::mutate(!!rlang::sym(sample_id) := as.character(!!rlang::sym(sample_id)))
  #- Add rownames
  rownames(design_mat_cln) <- design_mat_cln %>%
    pull(as.name(sample_id))
  #- Sample names
  cols_for_analysis <- design_mat_cln %>%
    pull(as.name(sample_id))

  #----- Counts Table
  # row_id %in% colnames(evidence_tbl_filt)
  # which(str_detect(setdiff(colnames(evidence_tbl_filt), row_id), group_pattern))
  evidence_tbl_filt <- vroom::vroom(counts_table_file, delim = "\t") %>%
    dplyr::select(one_of(c(row_id)), matches(group_pattern))
  #- Filter counts table for Sample Names
  evidence_tbl_col <- evidence_tbl_filt[, c(row_id, cols_for_analysis)]
  #- Save raw counts
  vroom::vroom_write(
    evidence_tbl_col,
    file.path(output_dir, "raw_counts.tsv")
  )
  writexl::write_xlsx(
    evidence_tbl_col,
    file.path(output_dir, "raw_counts.xlsx")
  )

  #----- Remove empty proteins without abundance data
  removeEmptyRows <- function(input_table, col_pattern, row_id) {
    #- temp_col_name <- "temp_uniprot_acc"
    temp_col_name <- paste0("temp_", as_name(enquo(row_id)))
    #- Add row numbers to column
    temp_input_table <- input_table %>%
      dplyr::mutate(!!rlang::sym(temp_col_name) := row_number())
    #- Row numbers that have atleast one count
    sites_to_accept <- temp_input_table %>%
      mutate(across(matches(col_pattern, perl = TRUE), \(x){
        (is.na(x) | x == 0)
      })) %>%
      dplyr::filter(!if_all(matches(col_pattern, perl = TRUE), \(x){
        x == TRUE
      })) %>%
      dplyr::select({{ temp_col_name }})
    #- Removing entries where all the "Reporter intensity corrected" rows are zero
    filtered_table <- temp_input_table %>%
      inner_join(sites_to_accept, by = temp_col_name) %>%
      dplyr::select(-temp_col_name)
    return(filtered_table)
  }
  #- Counts Table with atleast one count
  cln_dat_wide_unsorted <- removeEmptyRows(
    input_table = evidence_tbl_col,
    col_pattern = group_pattern,
    row_id = !!rlang::sym(row_id)
  )
  vroom::vroom_write(
    cln_dat_wide_unsorted,
    file.path(output_dir, "raw_counts_after_removing_empty_rows.tsv")
  )

  writexl::write_xlsx(
    cln_dat_wide_unsorted,
    file.path(output_dir, "raw_counts_after_removing_empty_rows.xlsx")
  )

  #----- Count the total number of missing values in total
  table_value <- cln_dat_wide_unsorted[, c(colnames(cln_dat_wide_unsorted)[1], cols_for_analysis)] %>%
    column_to_rownames(row_id) %>%
    log2() %>%
    data.matrix() %>%
    is.infinite() %>%
    table()
  #- Plot missing values
  plotNumMissingValues <- function(input_table) {
    plot_num_missing_values <- apply(
      data.matrix(log2(input_table)), 2,
      function(x) {
        length(which(!is.finite(x)))
      }
    ) |>
      t() |>
      t() |>
      set_colnames("No. of Missing Values") |>
      as.data.frame() |>
      rownames_to_column("Samples ID") |>
      ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90))
    plot_num_missing_values
  }
  plot_num_missing_values_before <- plotNumMissingValues(cln_dat_wide_unsorted[, cols_for_analysis])
  #- save plots
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("num_missing_values_before_filtering.", format_ext)
    )
    ggsave(
      filename = file_name,
      plot = plot_num_missing_values_before,
      limitsize = FALSE
    )
  }

  #----- Remove Rows with Missing Values per group
  #- For each experimental group, identify proteins that have more than accepted number of missing values per group
  removeRowsWithMissingValues <- function(input_table,
                                          cols,
                                          design_matrix,
                                          sample_id,
                                          row_id,
                                          group_column,
                                          max_num_samples_miss_per_group,
                                          abundance_threshold,
                                          temporary_abundance_column = "Abundance") {
    abundance_long <- input_table |>
      pivot_longer(
        cols = {{ cols }},
        names_to = as_name(enquo(sample_id)),
        values_to = temporary_abundance_column
      ) |>
      mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)) |>
      left_join(
        design_matrix |>
          mutate({{ sample_id }} := purrr::map_chr({{ sample_id }}, as.character)),
        by = as_name(enquo(sample_id))
      )

    count_missing_values_per_group <- abundance_long |>
      mutate(is_missing = ifelse(!is.na(!!sym(temporary_abundance_column)) &
        !!sym(temporary_abundance_column) > abundance_threshold, 0, 1)) |>
      group_by({{ row_id }}, {{ group_column }}) |>
      summarise(num_missing_values = sum(is_missing)) |>
      ungroup()

    remove_rows_temp <- count_missing_values_per_group |>
      dplyr::filter(max_num_samples_miss_per_group < num_missing_values) |>
      dplyr::select(-num_missing_values, -{{ group_column }}) |>
      distinct({{ row_id }})

    filtered_tbl <- input_table |>
      dplyr::anti_join(remove_rows_temp,
        by = as_name(enquo(row_id))
      )
    return(filtered_tbl)
  }
  #- Remove rows with missing values per group
  if (!is.na(max_num_samples_miss_per_group)) {
    cln_dat_wide_cleaned <- removeRowsWithMissingValues(
      cln_dat_wide_unsorted,
      matches(group_pattern),
      design_mat_cln |>
        dplyr::mutate(Sample_ID = as.character(sample_id)),
      !!rlang::sym(sample_id),
      !!rlang::sym(row_id),
      !!rlang::sym(group_id),
      max_num_samples_miss_per_group,
      abundance_threshold
    )
  } else {
    cln_dat_wide_cleaned <- cln_dat_wide_unsorted
  }
  #- Remove columns that dont match
  cln_dat_wide <- cln_dat_wide_cleaned[, c(colnames(cln_dat_wide_cleaned)[1], cols_for_analysis)]
  vroom::vroom_write(
    cln_dat_wide,
    file.path(
      output_dir,
      "raw_counts_after_removing_proteins_with_missing_values.tsv"
    )
  )
  writexl::write_xlsx(
    cln_dat_wide,
    file.path(
      output_dir,
      "raw_counts_after_removing_proteins_with_missing_values.xlsx"
    )
  )
  print(paste("Number of proteins before removing proteins with missing values:", nrow(evidence_tbl_col)))
  print(paste("Number of protein after proteins with missing values were removed:", nrow(cln_dat_wide)))
  print(paste("Number of proteins removed:", (nrow(evidence_tbl_col) - nrow(cln_dat_wide))))

  #----- Count the total number of missing values in total
  counts_filt <- cln_dat_wide |>
    column_to_rownames(row_id)
  #- plot missing values after filtering
  plot_num_missing_values <- plotNumMissingValues(counts_filt[, cols_for_analysis])
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("num_missing_values.", format_ext)
    )
    ggsave(
      filename = file_name,
      plot = plot_num_missing_values_before,
      limitsize = FALSE
    )
  }
  #- Count the number of values in table
  plotNumOfValues <- function(input_table) {
    plot_num_missing_values <- apply(
      data.matrix(log2(input_table)), 2,
      function(x) {
        length(which(!is.na(x)))
      }
    ) |>
      t() |>
      t() |>
      set_colnames("No. of Values") |>
      as.data.frame() |>
      rownames_to_column("Samples ID") |>
      ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90))
    plot_num_missing_values
  }
  plot_num_of_values <- plotNumOfValues(counts_filt[, cols_for_analysis])
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("num_of_values.", format_ext)
    )
    ggsave(filename = file_name, plot = plot_num_of_values, limitsize = FALSE)
  }
  #- Dealing with zero versus NA values
  counts_na <- counts_filt
  na_values_marker <- (counts_filt == 0)
  if (handling_zeros == "mark_as_na") {
    # Mark entries with zero values as NA, take log of values.
    counts_na[na_values_marker] <- NA
  } else if (handling_zeros == "add_pseudocount") {
    # Add pseudocount offset
    counts_na[na_values_marker] <- counts_na[na_values_marker] + .Machine$double.xmin
  }
  counts_na.log <- log2(counts_na)

  #----- Between Array Normalization
  counts_na.log.quant <- normalizeBetweenArrays(counts_na.log, method = normalization)
  vroom::vroom_write(
    as.data.frame(counts_na.log.quant) |>
      rownames_to_column(row_id),
    file.path(
      output_dir,
      "counts_after_normalization_before_imputation.tsv"
    )
  )
  writexl::write_xlsx(
    as.data.frame(counts_na.log.quant) |>
      rownames_to_column(row_id),
    file.path(
      output_dir,
      "counts_after_normalization_before_imputation.xlsx"
    )
  )

  #----- Calculate average value from technical replicates
  #- clean design matrix
  design_mat_updated <- design_mat_cln |>
    dplyr::filter(!is.na(!!rlang::sym(group_id)))
  #- average technical replicates
  counts_rnorm.log.for.imputation <- counts_na.log.quant
  averageValuesFromReplicates <- function(input_table,
                                          design_matrix,
                                          group_pattern,
                                          row_id,
                                          sample_id,
                                          average_replicates_id) {
    output_table <- input_table |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        names_to = sample_id,
        values_to = "value"
      ) |>
      left_join(design_matrix,
        by = sample_id
      ) |>
      group_by(!!rlang::sym(average_replicates_id), !!rlang::sym(row_id)) |>
      summarise(value = mean(value, na.rm = TRUE)) |>
      ungroup() |>
      pivot_wider(
        names_from = !!rlang::sym(average_replicates_id),
        values_from = "value"
      ) |>
      column_to_rownames(row_id) |>
      as.matrix()
    return(output_table)
  }

  if (!is.na(average_replicates_id)) {
    #- Average replicate
    counts_rnorm.log.for.imputation <- averageValuesFromReplicates(
      counts_na.log.quant,
      design_mat_cln,
      group_pattern,
      row_id,
      sample_id,
      average_replicates_id
    )
    #- Save output
    vroom::vroom_write(
      as.data.frame(counts_rnorm.log.for.imputation) |>
        rownames_to_column(row_id),
      file.path(
        output_dir,
        "counts_after_normalization_before_imputation_averaged.tsv"
      )
    )
    writexl::write_xlsx(
      as.data.frame(counts_rnorm.log.for.imputation) |>
        rownames_to_column(args$row_id),
      file.path(
        output_dir,
        "counts_after_normalization_before_imputation_averaged.xlsx"
      )
    )
    # Update Design
    design_mat_updated <- design_mat_cln |>
      mutate(!!rlang::sym(sample_id) := !!rlang::sym(average_replicates_id)) |>
      dplyr::select(one_of(sample_id, group_id)) |>
      dplyr::filter(!is.na(!!rlang::sym(group_id))) |>
      distinct()
    rownames(design_mat_updated) <- design_mat_updated[, sample_id]
    design_mat_updated <- design_mat_updated[colnames(counts_rnorm.log.for.imputation), ]
    vroom::vroom_write(
      design_mat_updated,
      file.path(
        output_dir,
        "design_matrix_prot_avg.tab"
      )
    )
  }

  #----- Heatmap
  plotNumOfValuesNoLog <- function(input_table) {
    plot_num_missing_values <- apply(
      data.matrix(input_table), 2,
      function(x) {
        length(which(!is.na(x)))
      }
    ) |>
      t() |>
      t() |>
      set_colnames("No. of Values") |>
      as.data.frame() |>
      rownames_to_column("Samples ID") |>
      ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90))
    plot_num_missing_values
  }
  plot_num_of_values_before_imputation <- plotNumOfValuesNoLog(counts_rnorm.log.for.imputation)
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("num_of_values_before_imputation.", format_ext)
    )
    ggsave(filename = file_name, plot = plot_num_of_values_before_imputation, limitsize = FALSE)
  }

  #----- Missing value imputation
  counts_rnorm.log.for.contrast <- NA
  impute.selectGrps.filtered <- NA
  count_num_replicates_per_protein <- NA

  if (imputation == TRUE) {
    imputation_groups <- data.frame(temp_id = colnames(counts_rnorm.log.for.imputation) |>
      str_split("_") |>
      purrr::map_chr(1)) |>
      left_join(design_mat_updated,
        by = c("temp_id" = sample_id)
      ) |>
      pull(!!sym(args$group_id))
    impute.selectGrps.filtered <- selectGrps(counts_rnorm.log.for.imputation,
      imputation_groups,
      percent = impute_min_percent,
      n = impute_min_num_of_groups
    )
    impute.scImpute.output <- scImpute(
      impute.selectGrps.filtered,
      impute_specific_percent,
      imputation_groups
    )[, colnames(impute.selectGrps.filtered)]
    imputed_values <- tImpute(impute.scImpute.output,
      assay = "imputed"
    )
    num_samples_per_group <- design_mat_updated |>
      group_by(!!sym(group_id)) |>
      summarise(num_samples_per_group = n()) |>
      ungroup()
    count_num_replicates_per_protein <- impute.selectGrps.filtered |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        values_to = "LogIntensity",
        names_to = sample_id
      ) |>
      left_join(design_mat_updated, by = sample_id) |>
      dplyr::filter(!is.na(LogIntensity)) |>
      group_by(!!sym(row_id), !!sym(group_id)) |>
      summarise(counts = n()) |>
      ungroup() |>
      left_join(num_samples_per_group,
        by = group_id
      ) |>
      dplyr::filter(counts / num_samples_per_group >= remove_imputed_perc) |>
      dplyr::select(-counts, -num_samples_per_group)
    vroom::vroom_write(
      as.data.frame(imputed_values) |>
        rownames_to_column(row_id),
      file.path(
        output_dir,
        "counts_after_normalization_and_imputation.tsv"
      )
    )
    writexl::write_xlsx(
      as.data.frame(imputed_values) |>
        rownames_to_column(row_id),
      file.path(
        output_dir,
        "counts_after_normalization_and_imputation.xlsx"
      )
    )
    if (length(which(is.na(imputed_values))) > 0 ||
      length(which(is.nan(imputed_values))) > 0) {
      stop("Imputation didn't work properly.")
    }
    counts_rnorm.log.for.contrast <- imputed_values
  } else {
    counts_rnorm.log.for.contrast <- counts_rnorm.log.for.imputation |>
      as.data.frame() |>
      as.matrix()
  }

  #----- Remove samples in which the design matrix has NA values
  list_of_sample_ids_to_use <- design_mat_updated |>
    dplyr::filter(!is.na(!!rlang::sym(group_id))) |>
    pull(!!rlang::sym(sample_id))

  counts_rnorm.log.for.contrast.na.rm <- counts_rnorm.log.for.contrast[, list_of_sample_ids_to_use]

  if (length(which(list_of_sample_ids_to_use != colnames(counts_rnorm.log.for.contrast.na.rm))) > 0) {
    stop("Number of samples in design matrix does not match number of samples in data matrix")
  }

  #----- Imputation
  if (imputation == TRUE & remove_imputed == TRUE) {
    imputed_values_remove_imputed <- counts_rnorm.log.for.contrast.na.rm |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        values_to = "LogIntensity",
        names_to = sample_id
      ) |>
      left_join(design_mat_updated,
        by = sample_id
      ) |>
      dplyr::filter(!is.na(LogIntensity)) |>
      dplyr::inner_join(count_num_replicates_per_protein,
        by = c(row_id, group_id)
      ) |>
      dplyr::select(-one_of(c(group_id))) |>
      arrange(row_id, sample_id) |>
      pivot_wider(
        id_cols = row_id,
        names_from = sample_id,
        values_from = "LogIntensity"
      )
    vroom::vroom_write(
      as.data.frame(imputed_values_remove_imputed),
      file.path(
        output_dir,
        "counts_after_normalization_and_remove_imputed.tsv"
      )
    )
    writexl::write_xlsx(
      as.data.frame(imputed_values_remove_imputed),
      file.path(
        output_dir,
        "counts_after_normalization_and_remove_imputed.xlsx"
      )
    )
  }

  #----- Run statistical tests without RUV
  ID <- rownames(counts_rnorm.log.for.contrast.na.rm)
  #- Assign experimental group list
  getTypeOfGrouping <- function(design_matrix,
                                group_id,
                                sample_id) {
    temp_type_of_grouping <- design_matrix |>
      dplyr::select(!!rlang::sym(group_id), !!rlang::sym(sample_id)) |>
      group_by(!!rlang::sym(group_id)) |>
      summarise(!!rlang::sym(sample_id) := list(!!rlang::sym(sample_id))) |>
      ungroup()
    type_of_grouping <- temp_type_of_grouping |> pull(!!rlang::sym(sample_id))
    names(type_of_grouping) <- temp_type_of_grouping |> pull(!!rlang::sym(group_id))
    return(type_of_grouping)
  }
  type_of_grouping <- getTypeOfGrouping(
    design_matrix = design_mat_cln,
    group_id = group_id,
    sample_id = sample_id
  )
  #- Run the linear model fitting and statistical tests for a set of contrasts,
  #- then adjust with Empirical Bayes function
  runTestsContrasts <- function(data,
                                contrast_strings,
                                design_matrix,
                                formula_string,
                                p_value_column = p.mod,
                                q_value_column = q.mod,
                                fdr_value_column = fdr.mod,
                                weights = NA,
                                treat_lfc_cutoff = NA,
                                eBayes_trend = FALSE,
                                eBayes_robust = FALSE) {
    ff <- as.formula(formula_string)
    mod_frame <- model.frame(ff, design_matrix)
    design_m <- model.matrix(ff, mod_frame)
    data_subset <- data[, rownames(design_m)]
    #- Make contrasts
    contr.matrix <- makeContrasts(
      contrasts = contrast_strings,
      levels = colnames(design_m)
    )
    #- Attach weights
    if (!is.na(weights)) {
      if (nrow(weights) == nrow(design_m)) {
        design_m <- cbind(design_m, weights)
      } else {
        stop("Stop: nrow(weights) should be equal to nrow(design_m)")
      }
    }
    #- Model
    fit <- lmFit(data_subset, design = design_m)
    cfit <- contrasts.fit(fit, contrasts = contr.matrix)
    eb.fit <- eBayes(cfit, trend = eBayes_trend, robust = eBayes_robust)
    #- Run treat over here
    t.fit <- NA
    result_tables <- NA
    #- assign log fold change threshold below which is scientifically not relevant
    if (!is.na(treat_lfc_cutoff)) {
      t.fit <- treat(eb.fit, lfc = as.double(treat_lfc_cutoff))
      result_tables <- purrr::map(contrast_strings, function(contrast) {
        de_tbl <- topTreat(t.fit, coef = contrast, n = Inf) |>
          mutate({{ q_value_column }} := qvalue(P.Value)$q) |>
          mutate({{ fdr_value_column }} := p.adjust(P.Value, method = "BH")) |>
          dplyr::rename({{ p_value_column }} := P.Value)
      })
    } else {
      t.fit <- eb.fit
      result_tables <- purrr::map(contrast_strings, function(contrast) {
        de_tbl <- topTable(t.fit, coef = contrast, n = Inf) |>
          mutate({{ q_value_column }} := qvalue(P.Value)$q) |>
          mutate({{ fdr_value_column }} := p.adjust(P.Value, method = "BH")) |>
          dplyr::rename({{ p_value_column }} := P.Value)
      })
    }
    names(result_tables) <- contrast_strings
    return(list(results = result_tables, fit.eb = t.fit))
  }
  list_rnorm.log.quant.ruv.r0 <- NA
  myRes_rnorm.log.quant <- NA
  list_rnorm.log.quant.ruv.r0 <- runTestsContrasts(counts_rnorm.log.for.contrast.na.rm,
    contrast_strings = contrasts_tbl[, 1][[1]],
    design_matrix = design_mat_updated,
    formula_string = formula_string,
    weights = NA,
    treat_lfc_cutoff = as.double(treat_lfc_cutoff),
    eBayes_trend = as.logical(eBayes_trend),
    eBayes_robust = as.logical(eBayes_robust)
  )
  myRes_rnorm.log.quant <- list_rnorm.log.quant.ruv.r0$results
  #- mean-variance relationship plots
  pdf(file.path(output_dir, "plotSA_before_ruvIII.pdf"))
  plotSA(list_rnorm.log.quant.ruv.r0$fit.eb)
  dev.off()
  png(file.path(output_dir, "plotSA_before_ruvIII.png"))
  plotSA(list_rnorm.log.quant.ruv.r0$fit.eb)
  dev.off()

  #----- Find the list of negative control genes using ANOVA
  #- Identify negative control proteins for use in removal of unwanted variation, using an ANOVA test
  getNegCtrlProtAnova <- function(data_matrix,
                                  design_matrix,
                                  group_column = "group",
                                  num_neg_ctrl = 500,
                                  q_val_thresh = 0.05) {
    #- Inspired by matANOVA function from PhosR package:
    #- http://www.bioconductor.org/packages/release/bioc/html/PhosR.html
    grps <- design_matrix[colnames(data_matrix), group_column]
    ps <- apply(data_matrix, 1, function(x) {
      if (length(unique(grps[!is.na(x)])) > 1) {
        summary(stats::aov(as.numeric(x) ~ grps))[[1]][["Pr(>F)"]][1]
      } else {
        return(NA_real_)
      }
    })
    ps[is.na(ps)] <- 1
    aov <- qvalue(ps)$qvalues # added lambda = 0 to prevent Error in pi0est(p, ...)
    filtered_list <- aov[aov > q_val_thresh]
    list_size <- ifelse(num_neg_ctrl > length(filtered_list), length(filtered_list), num_neg_ctrl)
    control_genes <- names(sort(filtered_list, decreasing = TRUE)[1:list_size])
    # nrow(data_matrix) - length(control_genes)
    control_genes_index <- rownames(data_matrix) %in% control_genes
    names(control_genes_index) <- rownames(data_matrix)
    return(control_genes_index)
  }
  control_genes_index <- getNegCtrlProtAnova(counts_rnorm.log.for.contrast.na.rm,
    design_matrix = design_mat_updated,
    group_column = group_id,
    num_neg_ctrl = num_neg_ctrl,
    q_val_thresh = control_genes_q_val_thresh
  )
  vroom::vroom_write(
    data.frame(
      temp_col = names(control_genes_index),
      is_control_genes = control_genes_index
    ) |>
      set_colnames(c(row_id, "is_control_genes")),
    file.path(
      output_dir,
      "ctrl_genes_list_ruv3.tsv"
    ),
    delim = "\t"
  )
  writexl::write_xlsx(
    data.frame(
      temp_col = names(control_genes_index),
      is_control_genes = control_genes_index
    ) |>
      set_colnames(c(row_id, "is_control_genes")),
    file.path(
      output_dir,
      "ctrl_genes_list_ruv3.xlsx"
    )
  )
  print(paste("Total Number of Genes: ", length(control_genes_index)))
  print(paste("Number of Control Genes: ", length(which(control_genes_index))))

  #----- Draw canonical correlation plot
  ruv_groups <- data.frame(temp_column = colnames(counts_rnorm.log.for.contrast.na.rm)) |>
    dplyr::rename(!!rlang::sym(sample_id) := "temp_column") |>
    left_join(
      design_mat_updated |>
        dplyr::mutate(!!rlang::sym(sample_id) := as.character(!!rlang::sym(sample_id))),
      by = sample_id
    )

  cancorplot_r1 <- ruv_cancorplot(
    Y = t(counts_rnorm.log.for.contrast.na.rm),
    X = ruv_groups |>
      dplyr::filter(!is.na(!!rlang::sym(group_id))) |>
      pull(!!rlang::sym(group_id)),
    ctl = control_genes_index
  )
  cancorplot_r1 <- cancorplot_r1 +
    geom_vline(
      xintercept = ruv_k,
      linetype = "dotted",
      color = "red"
    )

  for (format_ext in plots_format) {
    file_name <- file.path(output_dir, paste0("cancor_plot_round_1.", format_ext))
    ggsave(plot = cancorplot_r1, filename = file_name, limitsize = FALSE)
  }

  #----- Run RUVIII on the input counts table
  #- Converts a design matrix to a biological replicate matrix for use with ruvIII
  getRuvIIIReplicateMatrix <- function(design_matrix,
                                       sample_id_column,
                                       group_column,
                                       temp_column = is_replicate_temp) {
    ruvIII_replicates_matrix <- design_matrix |>
      dplyr::select({{ sample_id_column }}, {{ group_column }}) |>
      mutate({{ temp_column }} := 1) |>
      pivot_wider(
        id_cols = {{ sample_id_column }},
        names_from = {{ group_column }},
        values_from = {{ temp_column }},
        values_fill = 0
      ) |>
      column_to_rownames(as_name(enquo(sample_id_column))) |>
      as.matrix()
    ruvIII_replicates_matrix
  }

  ruvIII_replicates_matrix <- getRuvIIIReplicateMatrix(
    design_mat_updated,
    !!rlang::sym(sample_id),
    !!rlang::sym(replicate_group_id)
  )

  #- edited from missMethyl source code https://rdrr.io/bioc/missMethyl/src/R/RUVfunctions.R
  cmriRUVfit <- function(Y, X, ctl,
                         Z = 1,
                         k = NULL,
                         method = c("inv", "rinv", "ruv4", "ruv3", "ruv2"),
                         M = NULL, ...) {
    method <- match.arg(method)
    if ((method %in% c("ruv4", "ruv3", "ruv2")) & is.null(k)) {
      stop("'k' cannot be NULL if method is 'ruv4', 'ruv3' or 'ruv2'.")
    }
    if (mode(ctl) != "logical") {
      stop("'ctl' must be a logical vector.")
    }
    if (is.data.frame(Y)) {
      Y <- data.matrix(Y)
    }
    if (mode(Y) != "numeric") {
      stop("'Y' must be a numeric matrix.")
    }
    if (method == "ruv3" & is.null(M)) {
      stop("'M' cannot be NULL if method is 'ruv3'")
    }
    Y <- t(Y)
    fit <- switch(method,
      inv = ruv::RUVinv(Y = Y, X = X, ctl = ctl, Z = Z, ...),
      rinv = ruv::RUVrinv(Y = Y, X = X, ctl = ctl, Z = Z, k = k, ...),
      ruv4 = ruv::RUV4(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...),
      ruv2 = ruv::RUV2(Y = Y, X = X, ctl = ctl, k = k, Z = Z, ...),
      ruv3 = ruv::RUVIII(Y = Y, M = M, ctl = ctl, k = k, ...)
    )
    return(fit)
  }
  counts_rnorm.log.ruvIII_v1 <- cmriRUVfit(
    Y = counts_rnorm.log.for.contrast.na.rm,
    X = ruv_groups$group,
    ctl = control_genes_index,
    Z = 1,
    k = ruv_k,
    method = ruv_method,
    M = ruvIII_replicates_matrix
  ) |> t()

  vroom::vroom_write(
    counts_rnorm.log.ruvIII_v1 |>
      as.data.frame() |>
      rownames_to_column(row_id),
    file.path(
      output_dir,
      "normalized_counts_after_ruv.tsv"
    )
  )
  writexl::write_xlsx(
    counts_rnorm.log.ruvIII_v1 |>
      as.data.frame() |>
      rownames_to_column(row_id),
    file.path(
      output_dir,
      "normalized_counts_after_ruv.xlsx"
    )
  )

  #----- Imputation After RUV3
  if (imputation == TRUE & remove_imputed == TRUE) {
    imputed_ruv_remove_imputed <- counts_rnorm.log.ruvIII_v1 |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        values_to = "LogIntensity",
        names_to = sample_id
      ) |>
      left_join(design_mat_updated,
        by = sample_id
      ) |>
      dplyr::filter(!is.na(LogIntensity)) |>
      dplyr::inner_join(count_num_replicates_per_protein,
        by = c(row_id, group_id)
      ) |>
      dplyr::select(-one_of(c(group_id))) |>
      arrange(row_id, sample_id) |>
      pivot_wider(
        id_cols = row_id,
        names_from = sample_id,
        values_from = "LogIntensity"
      )
    vroom::vroom_write(
      as.data.frame(imputed_ruv_remove_imputed),
      file.path(
        output_dir,
        "normalized_counts_after_ruv_remove_imputed.tsv"
      )
    )
    writexl::write_xlsx(
      as.data.frame(imputed_ruv_remove_imputed),
      file.path(
        output_dir,
        "normalized_counts_after_ruv_remove_imputed.xlsx"
      )
    )
  }

  #----- Draw the RLE and PCA plots
  plot_width <- 15
  plot_height <- 14
  rle_pca_plots_arranged <- NA

  plotRle <- function(Y,
                      rowinfo = NULL,
                      probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                      ylim = c(-0.5, 0.5)) {
    #  checks = check.ggplot()
    # if (checks) {
    rle <- t(apply(t(Y) - apply(Y, 2, function(x) {
      median(x, na.rm = TRUE)
    }), 2, function(x) {
      quantile(x, probs = probs, na.rm = TRUE)
    }))
    colnames(rle) <- c("min", "lower", "middle", "upper", "max")
    df <- cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))
    if (!is.null(rowinfo)) {
      rowinfo <- data.frame(rowinfo = rowinfo)
      df_temp <- cbind(df, rowinfo)
      my.x.factor.levels <- df_temp |>
        arrange(rowinfo) |>
        distinct(rle.x.factor) |>
        pull(rle.x.factor)
      df <- df_temp |>
        mutate(rle.x.factor = factor(rle.x.factor,
          levels = my.x.factor.levels
        )) |>
        arrange(rowinfo)
    }

    rleplot <- ggplot(df, aes_string(x = "rle.x.factor")) +
      geom_boxplot(
        aes_string(
          lower = "lower", middle = "middle",
          upper = "upper", max = "max", min = "min"
        ),
        stat = "identity"
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90) # , axis.ticks.x = element_blank()
      ) +
      theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
      geom_hline(yintercept = 0) +
      coord_cartesian(ylim = ylim)
    if (!is.null(rowinfo)) {
      if (ncol(rowinfo) == 1) {
        rleplot <- rleplot + aes(fill = rowinfo) + labs(fill = "")
      }
    }
    return(rleplot)
  }

  plotPca <- function(data,
                      design_matrix,
                      sample_id_column = Sample_ID,
                      group_column = group,
                      label_column = {{ sample_id_column }},
                      title,
                      geom.text.size = 11,
                      ...) {
    pca.res <- pca(t(as.matrix(data)))
    proportion_explained <- pca.res$prop_expl_var
    temp_tbl <- pca.res$variates$X |>
      as.data.frame() |>
      rownames_to_column(as_name(enquo(sample_id_column))) |>
      left_join(design_matrix, by = as_name(enquo(sample_id_column)))
    unique_groups <- temp_tbl |>
      distinct({{ group_column }}) |>
      pull({{ group_column }})
    output <- temp_tbl |>
      ggplot(aes(PC1, PC2, col = {{ group_column }}, label = {{ label_column }})) +
      geom_point() +
      geom_text_repel(size = geom.text.size, show.legend = FALSE) +
      xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
      ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
      labs(title = title) +
      theme(legend.title = element_blank())
    output
  }

  rlePcaPlotList <- function(list_of_data_matrix,
                             list_of_design_matrix,
                             sample_id_column = Sample_ID,
                             group_column = group,
                             list_of_descriptions) {
    rle_list <- purrr::pmap(
      list(
        data_matrix = list_of_data_matrix,
        description = list_of_descriptions,
        design_matrix = list_of_design_matrix
      ),
      function(data_matrix,
               description,
               design_matrix) {
        plotRle(t(as.matrix(data_matrix)),
          rowinfo = design_matrix[
            colnames(data_matrix),
            as_name(enquo(group_column))
          ]
        ) +
          labs(title = description)
      }
    )
    pca_list <- purrr::pmap(
      list(
        data_matrix = list_of_data_matrix,
        description = list_of_descriptions,
        design_matrix = list_of_design_matrix
      ),
      function(data_matrix,
               description,
               design_matrix) {
        plotPca(data_matrix,
          design_matrix = design_matrix,
          sample_id_column = {{ sample_id_column }},
          group_column = {{ group_column }},
          title = description,
          cex = 7
        )
      }
    )
    list_of_plots <- c(rle_list, pca_list)
    rle_pca_plots_arranged <- ggarrange(
      plotlist = list_of_plots,
      nrow = 2,
      ncol = length(list_of_descriptions),
      common.legend = FALSE,
      legend = "bottom",
      widths = 10,
      heights = 10
    )
    rle_pca_plots_arranged
  }

  if (imputation == TRUE & remove_imputed == TRUE) {
    rle_pca_plots_arranged <-
      rlePcaPlotList(
        list_of_data_matrix = list(
          counts_rnorm.log.for.contrast.na.rm,
          imputed_ruv_remove_imputed |>
            column_to_rownames(row_id) |>
            as.matrix()
        ),
        list_of_design_matrix = list(
          design_mat_updated,
          design_mat_updated
        ),
        sample_id_column = !!rlang::sym(sample_id),
        group_column = !!rlang::sym(group_id),
        list_of_descriptions = list(
          "Before RUVIII",
          "After RUVIII"
        )
      )
  } else {
    rle_pca_plots_arranged <-
      rlePcaPlotList(
        list_of_data_matrix = list(
          counts_rnorm.log.for.contrast.na.rm,
          counts_rnorm.log.ruvIII_v1
        ),
        list_of_design_matrix = list(
          design_mat_updated,
          design_mat_updated
        ),
        sample_id_column = !!rlang::sym(sample_id),
        group_column = !!rlang::sym(group_id),
        list_of_descriptions = list(
          "Before RUVIII",
          "After RUVIII"
        )
      )
  }

  for (format_ext in plots_format) {
    file_name <- file.path(output_dir, paste0("rle_pca_plots.", format_ext))
    ggsave(
      plot = rle_pca_plots_arranged,
      filename = file_name,
      limitsize = FALSE,
      width = plot_width,
      height = plot_height
    )
  }

  #----- Compare the different experimental groups and obtain lists of differentially expressed proteins
  list_rnorm.log.quant.ruv.r1 <- NA
  myRes_rnorm.log.quant.ruv.r1 <- NA
  #- requires statmod library
  list_rnorm.log.quant.ruv.r1 <- runTestsContrasts(counts_rnorm.log.ruvIII_v1,
    contrast_strings = contrasts_tbl[, 1][[1]],
    design_matrix = design_mat_updated,
    formula_string = formula_string,
    weights = NA,
    treat_lfc_cutoff = as.double(treat_lfc_cutoff),
    eBayes_trend = as.logical(eBayes_trend),
    eBayes_robust = as.logical(eBayes_robust)
  )
  myRes_rnorm.log.quant.ruv.r1 <- list_rnorm.log.quant.ruv.r1$results
  # mean-variance relationship plot
  pdf(file.path(output_dir, "plotSA_after_ruvIII.pdf"))
  plotSA(list_rnorm.log.quant.ruv.r1$fit.eb)
  dev.off()
  png(file.path(output_dir, "plotSA_after_ruvIII.png"))
  plotSA(list_rnorm.log.quant.ruv.r1$fit.eb)
  dev.off()
  saveRDS(
    list_rnorm.log.quant.ruv.r1$fit.eb,
    file.path(output_dir, "fit.eb.RDS")
  )

  #----- Prepare data for drawing the volcano plots
  #- Format results table for use in volcano plots,
  #- counting number of significant proteins, p-values distribution histogram
  getSignificantData <- function(list_of_de_tables,
                                 list_of_descriptions,
                                 row_id = uniprot_acc,
                                 p_value_column = p.mod,
                                 q_value_column = q.mod,
                                 fdr_value_column = fdr.mod,
                                 log_q_value_column = lqm,
                                 log_fc_column = logFC,
                                 comparison_column = comparison,
                                 expression_column = expression,
                                 facet_column = analysis_type,
                                 q_val_thresh = 0.05) {
    get_row_binded_table <- function(de_table_list, description) {
      output <- purrr::map(
        de_table_list,
        function(tbl) {
          tbl |>
            rownames_to_column(as_name(enquo(row_id))) |>
            dplyr::select(
              {{ row_id }},
              {{ p_value_column }},
              {{ q_value_column }},
              {{ fdr_value_column }},
              {{ log_fc_column }}
            )
        }
      ) |>
        purrr::map2(names(de_table_list), \(.x, .y){
          .x |>
            mutate({{ comparison_column }} := .y)
        }) |>
        bind_rows() |>
        mutate({{ facet_column }} := description) |>
        separate({{ comparison_column }},
          sep = "=",
          into = c(
            as_name(enquo(comparison_column)),
            as_name(enquo(expression_column))
          )
        )
    }
    logfc_tbl_all <- purrr::map2(
      list_of_de_tables, list_of_descriptions,
      function(a, b) {
        get_row_binded_table(de_table_list = a, description = b)
      }
    ) |>
      bind_rows()
    selected_data <- logfc_tbl_all |>
      mutate({{ log_q_value_column }} := -log10(q.mod)) |>
      dplyr::select(
        {{ row_id }},
        {{ log_q_value_column }},
        {{ q_value_column }},
        {{ p_value_column }},
        {{ log_fc_column }},
        {{ comparison_column }},
        {{ expression_column }},
        {{ facet_column }}
      ) |>
      dplyr::mutate(colour = case_when(
        abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} >= q_val_thresh ~ "orange",
        abs({{ log_fc_column }}) >= 1 & {{ q_value_column }} < q_val_thresh ~ "purple",
        abs({{ log_fc_column }}) < 1 & {{ q_value_column }} < q_val_thresh ~ "blue",
        TRUE ~ "black"
      )) |>
      dplyr::mutate(colour = factor(colour,
        levels = c("black", "orange", "blue", "purple")
      ))
    selected_data
  }

  selected_data <- getSignificantData(
    list_of_de_tables = list(
      myRes_rnorm.log.quant,
      myRes_rnorm.log.quant.ruv.r1
    ),
    list_of_descriptions = list(
      "No RUV",
      "RUV applied"
    ),
    row_id = !!sym(row_id),
    p_value_column = p.mod,
    q_value_column = q.mod,
    fdr_value_column = fdr.mod,
    log_q_value_column = lqm,
    log_fc_column = logFC,
    comparison_column = comparison,
    expression_column = expression,
    facet_column = analysis_type,
    q_val_thresh = q_val_thresh
  ) |>
    dplyr::rename(log2FC = "logFC")

  #----- Remove rows with imputed values
  if (imputation & remove_imputed) {
    included_comparisons <- count_num_replicates_per_protein |>
      dplyr::inner_join(count_num_replicates_per_protein,
        by = c(row_id)
      ) |>
      dplyr::filter(group.x != group.y) |>
      dplyr::mutate(expression = paste0("group", group.x, "-group", group.y)) |>
      dplyr::select(-group.x, -group.y)
    join_condition <- rlang::set_names(
      c(row_id, "expression"),
      c(row_id, "expression")
    )
    selected_data <- selected_data |>
      mutate(q.mod.old = q.mod) |>
      dplyr::inner_join(included_comparisons,
        by = join_condition
      ) |>
      group_by(analysis_type, comparison, expression) |>
      nest() |>
      ungroup() |>
      mutate(data = purrr::map(data, function(x) {
        x |>
          mutate(q.mod = qvalue(p.mod)$qvalues)
      })) |>
      unnest(cols = c(data))
  }

  selected_data |>
    dplyr:::select(-colour, -lqm) |>
    vroom::vroom_write(file.path(
      output_dir,
      "lfc_qval_long.tsv"
    ))
  selected_data |>
    dplyr:::select(-colour, -lqm) |>
    writexl::write_xlsx(file.path(
      output_dir,
      "lfc_qval_long.xlsx"
    ))

  #----- Print the volcano plots
  #- Draw the volcano plot
  plotVolcano <- function(selected_data,
                          log_q_value_column = lqm,
                          log_fc_column = logFC,
                          q_val_thresh = 0.05,
                          formula_string = "analysis_type ~ comparison") {
    volplot_gg.all <- selected_data |>
      ggplot(aes(y = {{ log_q_value_column }}, x = {{ log_fc_column }})) +
      geom_point(aes(col = colour)) +
      scale_colour_manual(
        values = c(levels(selected_data$colour)),
        labels = c(
          paste0("Not significant, logFC > ", 1),
          paste0("Significant, logFC >= ", 1),
          paste0("Significant, logFC <", 1),
          "Not Significant"
        )
      ) +
      geom_vline(xintercept = 1, colour = "black", size = 0.2) +
      geom_vline(xintercept = -1, colour = "black", size = 0.2) +
      geom_hline(yintercept = -log10(q_val_thresh)) +
      theme_bw() +
      xlab("Log fold changes") +
      ylab("-log10 q-value") +
      theme(legend.position = "none")
    volplot_gg.plot <- volplot_gg.all
    if (!is.na(formula_string) | formula_string != "") {
      volplot_gg.plot <- volplot_gg.all +
        facet_grid(as.formula(formula_string),
          labeller = labeller(facet_category = label_wrap_gen(width = 10))
        )
    }
    volplot_gg.plot
  }
  volplot_gg.all <- plotVolcano(selected_data,
    log_q_value_column = lqm,
    log_fc_column = log2FC,
    q_val_thresh = q_val_thresh,
    formula_string = "analysis_type ~ comparison"
  )
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("volplot_gg_all.", format_ext)
    )
    ggsave(filename = file_name, plot = volplot_gg.all, width = 7.29, height = 6)
  }

  #----- Count the number of up or down significnat differentially expressed proteins
  #- Count the number of statistically significant differentially expressed proteins
  #- (according to user-defined threshold)
  countStatDeGenes <- function(data,
                               lfc_thresh = 0,
                               q_val_thresh = 0.05,
                               log_fc_column = log2FC,
                               q_value_column = q.mod) {
    # comparison <- as.data.frame(data) |>
    #   distinct(comparison) |>
    #   pull(comparison)
    selected_data <- data |>
      dplyr::mutate(status = case_when(
        {{ q_value_column }} >= q_val_thresh ~ "Not significant",
        {{ log_fc_column }} >= lfc_thresh &
          {{ q_value_column }} < q_val_thresh ~ "Significant and Up",
        {{ log_fc_column }} < lfc_thresh &
          {{ q_value_column }} < q_val_thresh ~ "Significant and Down",
        TRUE ~ "Not significant"
      ))
    counts <- selected_data |>
      group_by(status) |>
      summarise(counts = n()) |>
      ungroup()
    all_possible_status <- data.frame(status = c(
      "Not significant",
      "Significant and Up",
      "Significant and Down"
    ))
    results <- all_possible_status |>
      left_join(counts, by = c("status" = "status")) |>
      mutate(counts = ifelse(is.na(counts), 0, counts))
    return(results)
  }
  #- Format results table for use in volcano plots, counting number of significant proteins,
  #- p-values distribution histogram
  printCountDeGenesTable <- function(list_of_de_tables,
                                     list_of_descriptions,
                                     formula_string = "analysis_type ~ comparison",
                                     facet_column = analysis_type,
                                     comparison_column = comparison,
                                     expression_column = expression) {
    count_stat_de_genes_helper <- function(de_table, description) {
      purrr::map(de_table, \(x){
        countStatDeGenes(x,
          lfc_thresh = 0,
          q_val_thresh = 0.05,
          log_fc_column = logFC,
          q_value_column = q.mod
        )
      }) |>
        purrr::map2(names(de_table), \(.x, .y){
          .x |>
            mutate({{ comparison_column }} := .y)
        }) |>
        bind_rows() |>
        mutate({{ facet_column }} := description) |>
        separate({{ comparison_column }},
          sep = "=",
          into = c(
            as_name(enquo(comparison_column)),
            as_name(enquo(expression_column))
          )
        )
    }
    num_significant_de_genes_all <- purrr::map2(
      list_of_de_tables,
      list_of_descriptions,
      function(a, b) {
        count_stat_de_genes_helper(
          de_table = a,
          description = b
        )
      }
    ) |>
      bind_rows()
    num_sig_de_genes_barplot <- num_significant_de_genes_all |>
      dplyr::filter(status != "Not significant") |>
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = "identity", aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90))
    # print(head(num_sig_de_genes_barplot))
    if (!is.na(formula_string)) {
      num_sig_de_genes_barplot <- num_sig_de_genes_barplot +
        facet_grid(as.formula(formula_string))
    }
    return(list(plot = num_sig_de_genes_barplot, table = num_significant_de_genes_all))
  }
  num_sig_de_molecules <-
    printCountDeGenesTable(
      list_of_de_tables = list(
        myRes_rnorm.log.quant,
        myRes_rnorm.log.quant.ruv.r1
      ),
      list_of_descriptions = list(
        "No RUV",
        "RUV applied"
      ),
      formula_string = "analysis_type ~ comparison"
    )
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("num_sda_entities_barplot.", format_ext)
    )
    ggsave(
      filename = file_name,
      plot = num_sig_de_molecules$plot,
      height = 10,
      width = 7
    )
  }
  vroom::vroom_write(
    num_sig_de_molecules$table,
    file.path(
      output_dir,
      "num_significant_differentially_abundant_all.tab"
    )
  )
  writexl::write_xlsx(
    num_sig_de_molecules$table,
    file.path(
      output_dir,
      "num_significant_differentially_abundant_all.xlsx"
    )
  )

  #----- Print p-values distribution figure
  #- Draw the p-values distribution plot
  printPValuesDistribution <- function(selected_data,
                                       p_value_column = p.mod,
                                       formula_string = "is_ruv_applied ~ comparison") {
    breaks <- c(0, 0.001, 0.01, 0.05, seq(0.1, 1, by = 0.1))
    # after_stat(density)
    pvalhist <- ggplot(selected_data, aes({{ p_value_column }})) +
      theme(axis.title.y = element_blank()) +
      xlab("P-value") +
      geom_histogram(aes(y = after_stat(density)),
        breaks = breaks,
        position = "identity",
        color = "black"
      ) +
      geom_histogram(aes(y = after_stat(density)),
        breaks = breaks,
        position = "identity"
      )
    if (!is.na(formula_string)) {
      pvalhist <- pvalhist +
        facet_grid(as.formula(formula_string))
    }
    pvalhist
  }
  pvalhist <- printPValuesDistribution(selected_data,
    p_value_column = p.mod,
    formula_string = "analysis_type ~ comparison"
  )
  for (format_ext in plots_format) {
    file_name <- file.path(
      output_dir,
      paste0("p_values_distn.", format_ext)
    )
    ggsave(
      filename = file_name,
      plot = pvalhist,
      height = 10,
      width = 7
    )
  }

  #----- Create wide format output file
  norm_counts <- NA
  counts_table_to_use <- counts_rnorm.log.ruvIII_v1

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    set_colnames(paste0(colnames(counts_table_to_use), ".log2norm")) |>
    rownames_to_column(row_id)

  raw_counts <- counts_filt |>
    as.data.frame() |>
    set_colnames(paste0(colnames(counts_filt), ".raw")) |>
    rownames_to_column(row_id)

  de_proteins_wide <- selected_data |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(
      id_cols = c(!!sym(row_id)),
      names_from = c(comparison),
      names_sep = ":",
      values_from = c(log2FC, q.mod, p.mod)
    ) |>
    left_join(norm_counts, by = row_id) |>
    left_join(raw_counts, by = row_id) |>
    dplyr::arrange(across(matches("q.mod"))) |>
    distinct()

  vroom::vroom_write(
    de_proteins_wide,
    file.path(
      output_dir,
      paste0(file_prefix, "_wide.tsv")
    )
  )
  writexl::write_xlsx(
    de_proteins_wide,
    file.path(
      output_dir,
      paste0(file_prefix, "_wide.xlsx")
    )
  )

  #----- Create long format output file
  #- Create the de_protein_long and de_phos_long tables
  createDeResultsLongFormat <- function(lfc_qval_tbl,
                                        norm_counts_input_tbl,
                                        raw_counts_input_tbl,
                                        row_id,
                                        sample_id,
                                        group_id,
                                        group_pattern,
                                        design_matrix_norm,
                                        design_matrix_raw) {
    norm_counts <- norm_counts_input_tbl |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        names_to = sample_id,
        values_to = "log2norm"
      ) |>
      left_join(design_matrix_norm, by = sample_id) |>
      group_by(!!sym(row_id), !!sym(group_id)) |>
      arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
      mutate(replicate_number = paste0("log2norm.", row_number())) |>
      ungroup() |>
      pivot_wider(
        id_cols = c(!!sym(row_id), !!sym(group_id)),
        names_from = replicate_number,
        values_from = log2norm
      ) |>
      mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

    raw_counts <- raw_counts_input_tbl |>
      as.data.frame() |>
      rownames_to_column(row_id) |>
      pivot_longer(
        cols = matches(group_pattern),
        names_to = sample_id,
        values_to = "raw"
      ) |>
      left_join(design_matrix_raw, by = sample_id) |>
      group_by(!!sym(row_id), !!sym(group_id)) |>
      arrange(!!sym(row_id), !!sym(group_id), !!sym(sample_id)) |>
      mutate(replicate_number = paste0("raw.", row_number())) |>
      ungroup() |>
      pivot_wider(
        id_cols = c(!!sym(row_id), !!sym(group_id)),
        names_from = replicate_number,
        values_from = raw
      ) |>
      mutate({{ group_id }} := purrr::map_chr(!!sym(group_id), as.character))

    left_join_columns <- rlang::set_names(
      c(row_id, group_id),
      c(row_id, "left_group")
    )
    right_join_columns <- rlang::set_names(
      c(row_id, group_id),
      c(row_id, "right_group")
    )

    de_proteins_long <- lfc_qval_tbl |>
      dplyr::select(-lqm, -colour, -analysis_type) |>
      dplyr::mutate(expression = str_replace_all(expression, group_id, "")) |>
      separate(expression, sep = "-", into = c("left_group", "right_group")) |>
      left_join(norm_counts, by = left_join_columns) |>
      left_join(norm_counts,
        by = right_join_columns,
        suffix = c(".left", ".right")
      ) |>
      left_join(raw_counts, by = left_join_columns) |>
      left_join(raw_counts,
        by = right_join_columns,
        suffix = c(".left", ".right")
      ) |>
      arrange(comparison, q.mod, log2FC) |>
      distinct()
    de_proteins_long
  }

  counts_table_to_use <- counts_rnorm.log.ruvIII_v1

  de_proteins_long <- createDeResultsLongFormat(
    lfc_qval_tbl = selected_data |>
      dplyr::filter(analysis_type == "RUV applied"),
    norm_counts_input_tbl = counts_table_to_use,
    raw_counts_input_tbl = counts_filt,
    row_id = row_id,
    sample_id = sample_id,
    group_id = group_id,
    group_pattern = group_pattern,
    design_matrix_norm = design_mat_updated,
    design_matrix_raw = design_mat_updated
  )

  vroom::vroom_write(
    de_proteins_long,
    file.path(
      output_dir,
      paste0(file_prefix, "_long.tsv")
    )
  )
  writexl::write_xlsx(
    de_proteins_long,
    file.path(
      output_dir,
      paste0(file_prefix, "_long.xlsx")
    )
  )
}
