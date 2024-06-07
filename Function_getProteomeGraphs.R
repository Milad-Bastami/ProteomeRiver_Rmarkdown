#---------- Proteome Graphs
#----- Packages
# suppressPackageStartupMessages({
#   p_load(gridExtra)
#   p_load(Glimma)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir <- "results/proteomics/publication_graphs"
# #- Directory path for temporary files
# tmp_dir <- "results/proteomics/cache"
# #- Directory path for input data files
# input_dir <- "results/proteomics/de_proteins"
# #- Input file with the design matrix
# design_matrix_file <- "data/proteomics_design.tab"
# #- Input file with the design matrix, before technical replicates were averaged 
# #- and merged into one sample
# before_avg_design_matrix_file <- NA
# #- A string describing the sample ID. This must be a column that exists in the design matrix
# sample_id <- "Sample_ID"
# #- A string describing the replicate group ID. This must be a column that exists in the design matrix
# group_id <- "group"
# #- A string describing the row id
# row_id <- "uniprot_acc"
# #- q-value threshold below which a protein has statistically significant differetial expression
# q_val_thresh <- 0.05
# #- File path for the de_proteins_long_annot.tsv annotation file, an output from annot_prots_cmd.R
# de_proteins_long_file <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
# #- A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg]
# plots_format <- list("pdf","png")
# #- Print the top X gene names f
# top_x_gene_name <- 5
# #- Whether the analysis is proteomics or phosphoproteomics
# data_type <- "proteomics"
# #- The column which contains the log2FC value
# log2fc_column <- "log2FC"
# #- The column which contains the FDR value
# fdr_column <- "q.mod"

#----- Function
getProteomeGraphs <- function(
    tmp_dir,
    output_dir,
    input_dir,
    design_matrix_file,
    before_avg_design_matrix_file,
    sample_id,
    group_id,
    row_id,
    q_val_thresh,
    de_proteins_long_file,
    plots_format,
    top_x_gene_name,
    data_type,
    log2fc_column,
    fdr_column
) {
  
  #----- Create Results Directories
  createDirectoryIfNotExists <- function(file_path, mode = "0777") {
    #- Create directory recursively if it doesn't exist
    if (!file.exists(file_path)) {
      dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)
    }
  }
  createDirectoryIfNotExists(file_path = output_dir)
  createDirectoryIfNotExists(file_path = tmp_dir)
  
  #----- Input Data
  #- Before normalization
  counts_before_averaging <- NA
  if(!is.na(before_avg_design_matrix_file)) {
    counts_before_averaging_file <- file.path(input_dir, 
                                              "counts_after_normalization_before_imputation.tsv")
    counts_before_averaging <- vroom::vroom(counts_before_averaging_file)
  }
  
  #- After averaging (before RUVIII)
  counts_rnorm.log.quant <- NA
  file_a <- file.path(input_dir, "counts_after_normalization_and_imputation.tsv")
  file_b <- file.path(input_dir, "counts_after_normalization_before_imputation_averaged.tsv")
  file_c <- file.path(input_dir, "counts_after_normalization_before_imputation.tsv")
  counts_file <- file_c
  if(file.exists(file_a)) {
    counts_file <- file_a
  } else if (file.exists(file_b)) {
    counts_file <- file_b
  }
  counts_rnorm.log.quant <- vroom::vroom(counts_file)
  
  #- After averaging and RUVIII
  counts_rnorm.log.ruvIII <- NA
  if(file.exists(file.path(input_dir, 
                           "normalized_counts_after_ruv_remove_imputed.tsv"))) {
    counts_rnorm.log.ruvIII <- vroom::vroom(file.path(input_dir,  
                                                      "normalized_counts_after_ruv_remove_imputed.tsv"))
  } else if(file.exists(file.path(input_dir, 
                                  "normalized_counts_after_ruv.tsv"))) {
    counts_rnorm.log.ruvIII <- vroom::vroom(file.path(input_dir, 
                                                      "normalized_counts_after_ruv.tsv"))
  } else {
    print("Normalized counts after RUV file not found.")
  }
  
  #----- Read Design Matrix
  design_mat_cln <- vroom::vroom(design_matrix_file) %>%
    as.data.frame() %>%
    dplyr::mutate(!!rlang::sym(sample_id) := as.character(!!rlang::sym(sample_id)))
  rownames(design_mat_cln) <- design_mat_cln %>% pull(as.name(sample_id))
  #- design matrix when technical replicates has been averaged and merged into one
  before_avg_design_mat_cln <- NA
  if(!is.na(before_avg_design_matrix_file)) {
    before_avg_design_mat_cln <- vroom::vroom(before_avg_design_matrix_file) %>%
      as.data.frame() %>%
      dplyr::mutate(!!rlang::sym(sample_id) := as.character(!!rlang::sym(sample_id)))
    rownames(before_avg_design_mat_cln) <- before_avg_design_mat_cln %>% pull(as.name(sample_id))
  }
  
  #----- RLE plots
  plotRle <- function(Y, 
                      rowinfo = NULL, 
                      probs = c(0.05, 0.25, 0.5, 0.75, 0.95), 
                      ylim = c(-0.5, 0.5)){
    #  checks = check.ggplot()
    # if (checks) {
    rle = t(apply(t(Y) - apply(Y, 2, function(x){
      median(x, na.rm=TRUE)}), 2, function(x){
        quantile(x, probs = probs, na.rm=TRUE)}))
    colnames(rle) = c("min", "lower", "middle", "upper", "max")
    df = cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))
    if (!is.null(rowinfo)) {
      rowinfo = data.frame(rowinfo = rowinfo)
      df_temp = cbind(df, rowinfo)
      my.x.factor.levels <- df_temp |>
        arrange(rowinfo) |>
        distinct(rle.x.factor) |>
        pull(rle.x.factor)
      df <- df_temp |>
        mutate(rle.x.factor = factor(rle.x.factor,
                                     levels = my.x.factor.levels)) |>
        arrange(rowinfo)
    }
    
    rleplot = ggplot(df, aes_string(x = "rle.x.factor")) +
      geom_boxplot(aes_string(lower = "lower", middle = "middle",
                              upper = "upper", max = "max", min = "min"),
                   stat = "identity") +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90) #, axis.ticks.x = element_blank()
      ) +
      theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
      geom_hline(yintercept = 0) +
      coord_cartesian(ylim = ylim)
    if (!is.null(rowinfo))
      if (ncol(rowinfo) == 1)
        rleplot = rleplot + aes(fill = rowinfo) + labs(fill = "")
    return(rleplot)
  }
  
  counts_rnorm.log.quant_mat <- counts_rnorm.log.quant %>%
    as.data.frame() %>%
    column_to_rownames(row_id)
  
  counts_rnorm.log.ruvIII_mat <- counts_rnorm.log.ruvIII %>%
    as.data.frame() %>%
    column_to_rownames(row_id)
  #- Before RUV3 RLE plot
  before_RUVIII_rle <- plotRle(t(as.matrix(counts_rnorm.log.quant_mat)),
                               rowinfo = design_mat_cln[colnames(counts_rnorm.log.quant_mat),
                                                        group_id]) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")
  createDirectoryIfNotExists(file.path(output_dir, "RLE"))
  file_name_part <- file.path(output_dir, "RLE", "before_RUVIII_rle.")
  
  gg_save_logging <- function(input_plot, 
                              file_name_part, 
                              plots_format, 
                              width=7, 
                              height=7) {
    for(format_ext in plots_format) {
      file_name <- paste0(file_name_part, format_ext)
      ggsave(plot=input_plot, 
             filename = file_name, 
             width=width, 
             height=height)
    }
  }
  gg_save_logging(before_RUVIII_rle, file_name_part, plots_format)
  
  #- After RUV3 RLE plot
  after_RUVIII_rle <- plotRle(t(as.matrix(counts_rnorm.log.ruvIII_mat)),
                              rowinfo = design_mat_cln[colnames(counts_rnorm.log.ruvIII_mat),
                                                       group_id]) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")
  file_name_part <- file.path(output_dir, "RLE", "after_RUVIII_rle.")
  gg_save_logging(after_RUVIII_rle, file_name_part, plots_format)
  
  #- After RUV3 RLE plot After Design Matrix average
  if(!is.na(before_avg_design_matrix_file)) {
    counts_before_averaging_mat <- counts_before_averaging %>%
      as.data.frame() %>%
      column_to_rownames("uniprot_acc")
    counts_before_averaging_rle <- 
      plotRle(t(as.matrix(counts_before_averaging_mat)), 
              rowinfo = before_avg_design_mat_cln[colnames(counts_before_averaging_mat), group_id]) + 
      theme(axis.text.x = element_text(size = 13)) +
      theme(axis.text.y = element_text(size = 13)) +
      theme(axis.title.x = element_text(size = 12)) +
      theme(axis.title.y = element_text(size = 12)) +
      theme(plot.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 12)) +
      theme(legend.title = element_text(size = 12)) +
      xlab("Samples")
    file_name_part <- file.path(output_dir, "RLE", "counts_before_averaging_rle.")
    gg_save_logging(counts_before_averaging_rle, file_name_part, plots_format)
  }
  
  #----- Hierarchical clustering
  # before_hclust <- hclust(dist(t(as.matrix(counts_rnorm.log.quant_mat))))
  # plot(before_hclust)
  # 
  # after_hclust <- hclust(dist(t(as.matrix(counts_rnorm.log.ruvIII_mat))), method="ward.D2")
  # plot(after_hclust)
  
  #----- Volcano plots
  selected_data <- vroom::vroom(file.path(de_proteins_long_file)) %>% 
    mutate(lqm = -log10(!!sym(fdr_column))) %>%
    dplyr::mutate(label = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "Not sig., logFC >= 1",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC >= 1",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC < 1",
      TRUE ~ "Not sig.")) %>% 
    dplyr::mutate(colour = case_when(
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "orange",
      abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "purple",
      abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "blue",
      TRUE ~ "black")) %>%
    dplyr::mutate(colour = factor(colour, levels = c("black", "orange", "blue", "purple")))
  
  createDirectoryIfNotExists(file.path(output_dir, "Volcano_Plots"))
  
  #- Draw the volcano plot, used in publication graphs
  plotOneVolcano <- function(input_data, 
                             input_title,
                             log_q_value_column = lqm,
                             log_fc_column = logFC,
                             points_type_label = label,
                             points_color = colour,
                             q_val_thresh = 0.05) {
    colour_tbl <- input_data |>
      distinct({{points_type_label}}, {{points_color}})
    
    colour_map <- colour_tbl |>
      pull({{points_color}}) |>
      as.vector()
    
    names(colour_map) <- colour_tbl |>
      pull({{points_type_label}})
    
    avail_labels <- input_data |>
      distinct({{points_type_label}}) |>
      pull({{points_type_label}})
    
    avail_colours <- colour_map[avail_labels]
    
    volcano_plot <-  input_data |>
      ggplot(aes(y = {{log_q_value_column}},
                 x = {{log_fc_column}} )) +
      geom_point(aes(col = label)) +
      scale_colour_manual(values = avail_colours) +
      geom_vline(xintercept = 1, colour = "black", size = 0.2) +
      geom_vline(xintercept = -1, colour = "black", size = 0.2) +
      geom_hline(yintercept = -log10(q_val_thresh)) +
      theme_bw() +
      xlab(expression(Log[2](`fold-change`))) +
      ylab(expression(-log[10](`q-value`))) +
      labs(title = input_title)+  # Remove legend title
      theme(legend.title = element_blank()) +
      # theme(legend.position = "none")  +
      theme(axis.text.x = element_text(size = 13))   +
      theme(axis.text.y = element_text(size = 13))  +
      theme(axis.title.x = element_text(size = 12))  +
      theme(axis.title.y = element_text(size = 12))  +
      theme(plot.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 12)) # +
    # theme(legend.title = element_text(size = 12))
    volcano_plot
  }
  
  list_of_volcano_plots <- selected_data %>%
    group_by(comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate(title = paste( comparison)) %>%
    mutate(plot = purrr:::map2(data, title, \(x,y) { 
      plotOneVolcano(x, y, log_fc_column = !!sym(log2fc_column))}))
  
  # list_of_volcano_plots %>% pull(plot)
  purrr::walk2(list_of_volcano_plots %>% pull(title),
               list_of_volcano_plots %>% pull(plot),
               ~{file_name_part <- file.path(output_dir, "Volcano_Plots", paste0(.x, "."))
               gg_save_logging ( .y, file_name_part, plots_format)})
  
  ggsave(
    filename = file.path(output_dir, "Volcano_Plots", "list_of_volcano_plots.pdf" ),
    plot = marrangeGrob((list_of_volcano_plots %>% pull(plot)), nrow=1, ncol=1),
    width = 7, height = 7
  )
  
  #----- Significant genes Bar plot
  num_sig_de_molecules <- selected_data %>%
    dplyr::mutate(status = case_when(
      !!sym(fdr_column) >= 0.05 ~ "Not significant",
      !!sym(log2fc_column) >= 0 & !!sym(fdr_column) < 0.05 ~ "Significant and Up",
      !!sym(log2fc_column) < 0 & !!sym(fdr_column) < 0.05 ~ "Significant and Down",
      TRUE ~ "Not significant")) %>% 
    group_by(comparison, status) %>% 
    summarise(counts = n()) %>%
    ungroup()
  
  formula_string <- ". ~ comparison"
  
  num_sig_de_genes_barplot <- num_sig_de_molecules %>%
    dplyr::filter(status != "Not significant") %>%
    ggplot(aes(x = status, y = counts)) +
    geom_bar(stat = "identity") +
    geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
    theme(axis.text.x = element_text(angle = 90))  +
    facet_grid(as.formula(formula_string))
  
  createDirectoryIfNotExists(file.path(output_dir, "NumSigDeMolecules"))
  
  vroom::vroom_write(num_sig_de_molecules,
                     file.path(output_dir, "NumSigDeMolecules", "num_sig_de_molecules.tab" ) )
  
  num_of_comparison <- num_sig_de_molecules |>
    distinct(comparison) |>
    nrow()
  
  ggsave(filename = file.path(output_dir, "NumSigDeMolecules", "num_sig_de_molecules.png" ),
         plot = num_sig_de_genes_barplot,
         height = 10,
         width = (num_of_comparison + 2) *7/6 )
  
  #----- PCA plots
  plotPca <- function(data,
                      design_matrix,
                      sample_id_column = Sample_ID,
                      group_column = group,
                      label_column = {{sample_id_column}},
                      title, 
                      geom.text.size=11,
                      ...){
    pca.res <- pca(t(as.matrix(data)))
    proportion_explained <- pca.res$prop_expl_var
    temp_tbl <- pca.res$variates$X |>
      as.data.frame() |>
      rownames_to_column(as_name(enquo(sample_id_column))) |>
      left_join(design_matrix, by = as_name(enquo(sample_id_column)))
    unique_groups <- temp_tbl |> distinct( {{group_column}}) |> pull( {{group_column}})
    output <- temp_tbl |>
      ggplot(aes(PC1, PC2, col = {{group_column}}, label = {{label_column}})) +
      geom_point() +
      geom_text_repel(size  = geom.text.size, show.legend=FALSE) +
      xlab( paste( "PC1 (", round(proportion_explained$X[["PC1"]]*100, 0),"%)", sep="")) +
      ylab( paste( "PC2 (", round(proportion_explained$X[["PC2"]]*100, 0),"%)", sep="")) +
      labs(title = title) +
      theme(legend.title = element_blank())
    output
  }
  
  createDirectoryIfNotExists(file.path(output_dir, "PCA"))
  #- PCA before RUV3
  counts_before_averaging_mat <- NA
  if(!is.na(before_avg_design_matrix_file)) {
    counts_before_averaging_mat <- counts_before_averaging  %>%
      as.data.frame() %>%
      column_to_rownames(row_id)
  }
  
  before_ruvIII_pca <- plotPca(counts_rnorm.log.quant_mat,
                               design_matrix = design_mat_cln,
                               sample_id_column = !!rlang::sym(sample_id),
                               group_column = !!rlang::sym(group_id),
                               title = "Before RUVIII",
                               geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "before_ruvIII_pca.")
  gg_save_logging(before_ruvIII_pca, file_name_part, plots_format)
  
  #- PCA after RUV3
  after_ruvIII_pca <- plotPca(counts_rnorm.log.ruvIII_mat,
                              design_matrix = design_mat_cln,
                              sample_id_column = !!rlang::sym(sample_id),
                              group_column = !!rlang::sym(group_id),
                              title = "After RUVIII", 
                              geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "after_ruvIII_pca.")
  gg_save_logging(after_ruvIII_pca, file_name_part, plots_format)
  
  #- PCA before RUV3 label group
  before_ruvIII_pca <- plotPca(counts_rnorm.log.quant_mat,
                               design_matrix = design_mat_cln,
                               sample_id_column = !!rlang::sym(sample_id),
                               group_column = !!rlang::sym(group_id),
                               label_column = !!rlang::sym(group_id),
                               title = "Before RUVIII",
                               geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "before_ruvIII_pca_group_labels.")
  gg_save_logging(before_ruvIII_pca, file_name_part, plots_format)
  
  #- PCA after RUV3 label group
  after_ruvIII_pca <- plotPca(counts_rnorm.log.ruvIII_mat,
                              design_matrix = design_mat_cln,
                              sample_id_column = !!rlang::sym(sample_id),
                              group_column = !!rlang::sym(group_id),
                              label_column = !!rlang::sym(group_id),
                              title = "After RUVIII", 
                              geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "after_ruvIII_pca_group_labels.")
  gg_save_logging (after_ruvIII_pca, file_name_part, plots_format)
  
  #- PCA before RUV3 no labels
  before_ruvIII_pca_no_labels <- plotPca(counts_rnorm.log.quant_mat,
                                         design_matrix = design_mat_cln %>% 
                                           mutate( sample_labels = ""),
                                         sample_id_column = !!rlang::sym(sample_id),
                                         group_column = !!rlang::sym(group_id),
                                         label_column = sample_labels,
                                         title = "Before RUVIII",
                                         geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "before_ruvIII_no_sample_labels.")
  gg_save_logging(before_ruvIII_pca_no_labels, file_name_part, plots_format)
  
  #- PCA after RUV3 no labels
  after_ruvIII_pca_no_labels <- plotPca(counts_rnorm.log.ruvIII_mat,
                                        design_matrix = design_mat_cln %>%
                                          mutate( sample_labels = ""),
                                        sample_id_column = !!rlang::sym(sample_id),
                                        group_column = !!rlang::sym(group_id),
                                        label_column = sample_labels,
                                        title = "After RUVIII", 
                                        geom.text.size = 7) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))
  
  file_name_part <- file.path(output_dir, "PCA", "after_ruvIII_no_sample_labels.")
  gg_save_logging(after_ruvIII_pca_no_labels, file_name_part, plots_format)
  
  #- PCA before sample averages
  if(!is.na(before_avg_design_matrix_file)) {
    counts_before_averaging_pca <- plotPca(counts_before_averaging_mat,
                                           design_matrix = before_avg_design_mat_cln,
                                           sample_id_column = !!rlang::sym(sample_id),
                                           group_column = !!rlang::sym(group_id),
                                           title = "After RUVIII", 
                                           geom.text.size = 7) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12)) +
      theme(axis.text.y = element_text(size = 12)) +
      theme(axis.title.x = element_text(size = 12)) +
      theme(axis.title.y = element_text(size = 12)) +
      theme(plot.title = element_text(size = 12)) +
      theme(legend.text = element_text(size = 12)) +
      theme(legend.title = element_text(size = 12))
    
    file_name_part <- file.path( args$output_dir, "PCA", "counts_before_averaging_pca.")
    gg_save_logging ( counts_before_averaging_pca, file_name_part, args$plots_format)
  }
  
  #----- Interactive Volcano Plot
  #- Interactive Volcano Plot proteomics
  getGlimmaVolcanoProteomics <- function(r_obj, 
                                         coef, 
                                         volcano_plot_tab, 
                                         uniprot_column = best_uniprot_acc, 
                                         gene_name_column = gene_name, 
                                         display_columns = c("PROTEIN_NAMES"), 
                                         output_dir) {
    if(coef <= ncol(r_obj$coefficients)) {
      best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:") |>
        purrr::map_chr(1)
      volcano_plot_tab_cln <- volcano_plot_tab |>
        dplyr::distinct({{uniprot_column}}, 
                        {{gene_name_column}}, 
                        pick(one_of(display_columns))) |>
        dplyr::rename(best_uniprot_acc = {{uniprot_column}}, 
                      gene_name = {{gene_name_column}})
      anno_tbl <- data.frame(uniprot_acc = rownames(r_obj@.Data[[1]]), 
                             best_uniprot_acc = best_uniprot_acc) |> 
        left_join(volcano_plot_tab_cln, 
                  by = c("best_uniprot_acc")) |>
        mutate(gene_name = case_when(is.na(gene_name) ~ best_uniprot_acc, 
                                     TRUE ~ gene_name) )
      gene_names <- anno_tbl |>
        pull(gene_name)
      rownames(r_obj@.Data[[1]]) <- gene_names
        ## fix html filename issue (remove '=')
        html_filename <- 
        strsplit(colnames(r_obj$coefficients)[coef], split = '=', fixed = T)[[1]][1] |>
        paste0('.html')
      htmlwidgets::saveWidget(widget = glimmaVolcano(r_obj, 
                                                     coef = coef, 
                                                     anno = anno_tbl, 
                                                     display.columns = display_columns), 
                              file = file.path(output_dir,  html_filename),  
                              selfcontained = TRUE)
    }
  }
  
  if(data_type  == "proteomics" && 
     file.exists(file.path(input_dir, "fit.eb.RDS")) && 
     file.exists(de_proteins_long_file)) {
    
    volcano_plot_tab <- vroom::vroom(de_proteins_long_file)  %>%
      mutate(lqm = -log10(!!sym(fdr_column))) |>
      dplyr::mutate(label = case_when(
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "Not sig., logFC >= 1", 
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC >= 1", 
        abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC < 1", 
        TRUE ~ "Not sig.")) |>
      dplyr::mutate(colour = case_when(
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "orange", 
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "purple", 
        abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "blue", 
        TRUE ~ "black")) |>
      dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:") |> 
                      purrr::map_chr(1)) |>
      dplyr::mutate(best_uniprot_acc = str_split(!!sym(row_id), ":") |> 
                      purrr::map_chr(1)) |>
      dplyr::mutate(analysis_type = comparison) |>
      dplyr::rename(PROTEIN_NAMES = "PROTEIN-NAMES") |>
      dplyr::select(best_uniprot_acc, lqm, !!sym(fdr_column), p.mod, !!sym(log2fc_column), 
                    comparison, label, colour, gene_name, `PROTEIN_NAMES`) |>
      dplyr::mutate(my_alpha = case_when(gene_name != "" ~ 1, TRUE ~ 0.5))
    
    r_obj <- readRDS(file.path(input_dir, "fit.eb.RDS"))
    
    output_dir1 <- file.path(output_dir, "Interactive_Volcano_Plots")
    createDirectoryIfNotExists(output_dir1)
    
    purrr::walk(seq_len(ncol(r_obj$coefficients)), 
                \(coef) { 
                  getGlimmaVolcanoProteomics(r_obj, 
                                             coef = coef, 
                                             volcano_plot_tab = volcano_plot_tab, 
                                             uniprot_column = best_uniprot_acc, 
                                             gene_name_column = gene_name, 
                                             display_columns = c("best_uniprot_acc", "PROTEIN_NAMES"),  
                                             output_dir = output_dir1)
                })
  }
  
  #- Interactive Volcano Plot phosphoproteomics
  getGlimmaVolcanoPhosphoproteomics <- function(r_obj, 
                                                coef, 
                                                volcano_plot_tab, 
                                                sites_id_column = sites_id, 
                                                sites_id_display_column = sites_id_short, 
                                                display_columns = c("sequence", "PROTEIN_NAMES"), 
                                                output_dir) {
    if(coef <= ncol(r_obj$coefficients)) {
      volcano_plot_tab_cln <- volcano_plot_tab |>
        dplyr::distinct({{sites_id_column}}, 
                        {{sites_id_display_column}}, 
                        pick(one_of(display_columns)))
      anno_tbl <- data.frame(sites_id = rownames(r_obj@.Data[[1]])) |>
        left_join(volcano_plot_tab_cln, 
                  by = join_by(sites_id == {{sites_id_column}}))
      sites_id_short_list <- anno_tbl |> 
        pull(sites_id_short)
      rownames( r_obj@.Data[[1]] ) <- sites_id_short_list
        # fix interactive plot htm file name issue (remove '=')
       html_filename <- strsplit(colnames(r_obj$coefficients)[coef], split = '=', fixed = T)[[1]][1] |>
        paste0('.html')
      htmlwidgets::saveWidget(widget = glimmaVolcano(r_obj, 
                                                     coef=coef, 
                                                     anno=anno_tbl, 
                                                     display.columns=display_columns), 
                              file = file.path(output_dir, html_filename),  
                              selfcontained = TRUE)
    }
  }
  
  if(data_type == "phosphoproteomics" 
     && file.exists(file.path(input_dir, "fit.eb.RDS"))
     && file.exists(de_proteins_long_file)) {
    
    merge_residue_position_lists <- function(residue, position) {
      residues_list <- str_split(residue, ";")[[1]]
      positions_list <- str_split(position, ";")[[1]]
      
      if(length(residues_list) != length(positions_list)) {
        stop(paste("Length not equal", residue, position))
      }
      
      purrr::map2_chr(residues_list, 
                      positions_list, 
                      function(r, p){ paste0(r, p) }) |>
        paste(collapse=";")
    }
    
    clean_first_positiion <- function(position) {
      first_position <- str_split(position, "\\|") |>
        purrr::map_chr(1) |>
        str_replace("\\(", "") |>
        str_replace("\\)", "")
      first_position
    }
    
    de_proteins_long_tbl <- vroom::vroom(de_proteins_long_file)
    
    volcano_plot_colour_points <- de_proteins_long_tbl %>%
      mutate(lqm = -log10(!!sym(fdr_column))) |>
      dplyr::mutate(label = case_when(
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "Not sig., logFC >= 1", 
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC >= 1", 
        abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "Sig., logFC < 1", 
        TRUE ~ "Not sig.")) |>
      dplyr::mutate(colour = case_when(
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= q_val_thresh ~ "orange", 
        abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < q_val_thresh ~ "purple", 
        abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < q_val_thresh ~ "blue", 
        TRUE ~ "black")) |>
      dplyr::mutate(analysis_type = comparison) |>
      dplyr::mutate(gene_name = str_split(UNIPROT_GENENAME, " |:") |> 
                      purrr::map_chr(1)) |> 
      dplyr::mutate(best_uniprot_acc = str_split(!!sym(row_id), ":" ) |> 
                      purrr::map_chr(1)) |>
      mutate(first_position = purrr::map_chr(position, clean_first_positiion)) |>
      mutate(merged_sites_residues = purrr::map2_chr(residue,  
                                                     first_position, 
                                                     \(residue, position) {
                                                       merge_residue_position_lists(residue, position)})) |>
      mutate(sites_id_short = paste0(gene_name, ":", merged_sites_residues, ":", best_uniprot_acc)) |>
      relocate(sites_id_short, .after="sites_id")
    
    volcano_plot_tab <- volcano_plot_colour_points |>
      dplyr::rename(PROTEIN_NAMES = "PROTEIN-NAMES") |>
      dplyr::select(sites_id, sites_id_short, best_uniprot_acc, lqm, !!sym(fdr_column), 
                    !!sym(log2fc_column), comparison, label, colour,  gene_name, 
                    sequence, `PROTEIN_NAMES`) |>
      dplyr::mutate(my_alpha = case_when(gene_name !=  "" ~ 1, TRUE ~ 0.5))
    
    r_obj <- readRDS( file.path(input_dir, "fit.eb.RDS"))
    
    output_dir1 <- file.path(output_dir, "Interactive_Volcano_Plots")
    createDirectoryIfNotExists(output_dir1)
    
    purrr::walk(seq_len(ncol(r_obj$coefficients)), 
                \(coef) { 
                  getGlimmaVolcanoPhosphoproteomics(r_obj, 
                                                    coef = coef, 
                                                    volcano_plot_tab = volcano_plot_tab, 
                                                    sites_id_column = sites_id, 
                                                    sites_id_display_column = sites_id_short, 
                                                    display_columns = c("sequence", "PROTEIN_NAMES"),  
                                                    output_dir = output_dir1)})
  }
}
