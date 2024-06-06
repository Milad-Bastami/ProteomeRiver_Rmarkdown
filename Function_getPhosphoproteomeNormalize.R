#---------- Phosphoproteome Normalize
#----- Packages
# suppressPackageStartupMessages({
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/norm_phos_by_prot_abundance"
# #- Directory path for temporary files
# tmp_dir <- "results/phosphoproteomics/cache"
# #- File with table of differentially expressed proteins
# proteins_file <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
# #- File with table of differentially abundant phosphosites
# phospho_file <- "results/phosphoproteomics/annot_phos/de_phos_long_annot.tsv"
# #- A string choosing which method to identify best host-protein from each protein group 
# #- [protein_rank, best_p_value, best_log_fc]
# protein_group_filter <- "protein_rank"
# #- A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg]
# plots_format <- list("pdf","png")

#----- Function
getPhosphoproteomeNormalize <- function(
    tmp_dir,
    output_dir,
    proteins_file,
    phospho_file,
    protein_group_filter,
    plots_format
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
  
  #----- Read phosphopeptides abundance data
  phospho_tbl_orig <- vroom::vroom(phospho_file)
  
  list_of_phospho_columns <- c("sites_id", "uniprot_acc", "gene_name", "position", 
                               "residue", "sequence", "q.mod", "fdr.mod", 
                               "p.mod", "log2FC", "comparison", "maxquant_row_ids")
  
  if(!"fdr.mod" %in% colnames(phospho_tbl_orig)) {
    list_of_phospho_columns <- list_of_phospho_columns[list_of_phospho_columns != "fdr.mod"]
  }
  
  phospho_tbl <- phospho_tbl_orig %>%
    dplyr::select(one_of(list_of_phospho_columns)) %>%
    dplyr::mutate(phos_row_ids = sites_id)
  
  phospho_cln <- phospho_tbl %>%
    separate_rows(uniprot_acc, position, residue, sequence, sep=":") %>%
    mutate(phospho_row_rank = row_number())
  
  #----- Read proteomics abundance data
  proteins_tbl_orig <-  vroom::vroom(proteins_file, delim = "\t")
  
  list_of_prot_columns <- c("uniprot_acc", "q.mod", "fdr.mod", "p.mod", 
                            "log2FC", "comparison", "maxquant_row_id")
  
  if(!"fdr.mod" %in% colnames(proteins_tbl_orig)) {
    list_of_prot_columns <- list_of_prot_columns[list_of_prot_columns != "fdr.mod"]
  }
  
  proteins_tbl <- proteins_tbl_orig %>%
    dplyr::select(one_of(list_of_prot_columns) ) %>%
    dplyr::rename(prot_maxquant_row_ids = "maxquant_row_id")
  
  proteins_cln <- proteins_tbl %>%
    separate_rows(uniprot_acc, sep=":") %>%
    mutate(protein_row_rank = row_number())
  
  #----- Uniprot List
  proteins_uniprot_list <- proteins_cln %>% 
    distinct(uniprot_acc) %>% 
    pull(uniprot_acc)
  
  phospho_uniprot_list <- phospho_cln %>% 
    distinct(uniprot_acc) %>% 
    pull(uniprot_acc)
  
  prot_phos_uniprot_list <- intersect(proteins_uniprot_list, 
                                      phospho_uniprot_list)
  
  print(proteins_uniprot_list %>% length())
  print(phospho_uniprot_list %>% length())
  print(prot_phos_uniprot_list %>% length())
  
  #----- Normalisation of the phosphopeptide abundance with the protein abundance
  #- A phosphopeptide can be matched to multiple proteins and each of those protein 
  #- may have their own abundance level
  #- Here I find distinct protein row ID and phosphosites row ID pairs
  join_protein_phosopho_keys <- proteins_cln %>%
    inner_join(phospho_cln,
               by = c("uniprot_acc" = "uniprot_acc",
                      "comparison" = "comparison"), 
               suffix=c(".prot", ".phos")) %>%
    distinct(comparison, prot_maxquant_row_ids, phos_row_ids)
  
  #- Function to find intersection between the uniprot accession of the phosphopeptide 
  #- and the uniprot accession of the protein used to do fold-change normalization
  intersectTwoUniprotList <- function(x, y) {  
    paste(intersect(str_split(x, ":")[[1]], str_split(y, ":")[[1]]), collapse=":") }
  
  #- Which array position in the list of proteins were the selected proteins
  getArrayPositionSelected <- function(x, y){ 
    which(str_split(x, ":")[[1]] %in% str_split(y, ":")[[1]]) }
  
  #- Once I have identified which protein was used for fold-change normalization, 
  #- I have pick out the position, residue, and peptide sequence from a list that 
  #- belongs to that protein. This is the function to do that array subsetting
  subsetPhosphositeDetails <- function(x,y) { 
    paste(str_split(x, ":")[[1]][ y], collapse=":") }
  
  #- I then have to find the unique protein row ID and phosphosites row ID pairs among the orignal
  #- protein groups table and phosphopeptide table, where there could be multiple possible host-proteins per phosphopeptide
  basic_data_shared <- join_protein_phosopho_keys %>%
    left_join(phospho_tbl, 
              by = c("phos_row_ids", "comparison")) %>%
    left_join(proteins_tbl, 
              by = c("prot_maxquant_row_ids", "comparison" ), 
              suffix = c(".phos", ".prot")) %>%
    mutate(uniprot_acc = purrr::map2_chr(uniprot_acc.phos, 
                                         uniprot_acc.prot, 
                                         intersectTwoUniprotList)) %>%
    mutate(array_pos = purrr::map2(uniprot_acc.phos, 
                                   uniprot_acc, 
                                   getArrayPositionSelected)) %>%
    mutate(gene_name = purrr::map2_chr(gene_name, 
                                       array_pos, 
                                       subsetPhosphositeDetails)) %>%
    mutate(position = purrr::map2_chr(position, 
                                      array_pos, 
                                      subsetPhosphositeDetails)) %>%
    mutate(residue = purrr::map2_chr(residue, 
                                     array_pos, 
                                     subsetPhosphositeDetails)) %>%
    #mutate(residue = purrr::map2_chr(residue, sequence, subsetPhosphositeDetails)) %>%
    dplyr::select(-uniprot_acc.phos, -uniprot_acc.prot) %>%
    distinct() %>%
    ## Calculate the new log fold-change and use the Fisher's method to calculate the updated p-value
    mutate(norm_phos_logFC = log2FC.phos - log2FC.prot) %>%
    mutate(adj_qmod.prot = ifelse(sign(log2FC.phos) == sign(log2FC.prot), 1-q.mod.prot, q.mod.prot)) %>%
    mutate(combined_q_mod = 1-pchisq(-2*( log(q.mod.phos) + log(adj_qmod.prot)), 2*2 )) %>%
    dplyr::mutate(status  = "Phos_and_Prot")
  
  list_of_data_shared_columns <- c("comparison", "norm_phos_logFC", "combined_q_mod", "combined_fdr_mod", 
                                   "sites_id", "uniprot_acc", "position", "residue", "sequence",
                                   "log2FC.phos", "q.mod.phos", "log2FC.prot", "q.mod.prot",
                                   "status", "phos_row_ids", "prot_maxquant_row_ids")
  
  if("fdr.mod" %in% colnames(proteins_cln) &
     "fdr.mod" %in% colnames(phospho_cln)) {
    basic_data_shared <- basic_data_shared %>%
      mutate(adj_fdrmod.prot = ifelse(sign(log2FC.phos) == sign(log2FC.prot), 1-fdr.mod.prot, fdr.mod.prot)) %>%
      mutate(combined_fdr_mod = 1-pchisq(-2*(log(fdr.mod.phos) + log(adj_fdrmod.prot)), 2*2 )) %>%
      dplyr::select(-adj_fdrmod.prot)
    list_of_data_shared_columns <- c("comparison",  "norm_phos_logFC", "combined_q_mod", "combined_fdr_mod", 
                                     "sites_id", "uniprot_acc", "position", "residue", "sequence",
                                     "log2FC.phos", "q.mod.phos", "fdr.mod.phos",
                                     "log2FC.prot", "q.mod.prot", "fdr.mod.prot",
                                     "status", "phos_row_ids", "prot_maxquant_row_ids")
  }
  
  if(!"combined_fdr_mod" %in% colnames(basic_data_shared)) { 
    list_of_data_shared_columns <- 
      list_of_data_shared_columns[!list_of_data_shared_columns %in% 
                                    c("combined_fdr_mod", "fdr.mod.phos", "fdr.mod.prot")]
  }
  
  basic_data_shared <- basic_data_shared %>%
    dplyr::select(one_of(list_of_data_shared_columns)) %>%
    arrange(comparison, combined_q_mod, norm_phos_logFC) %>%
    distinct() %>%
    as.data.frame
  
  #----- norm_phosphosite_lfc_minus_protein_lfc_with_protein_repeats
  #- Find all the phosphopeptides which did not have a matching protein for normalization
  #- Then carry over the log fold-change and p-value without changing it
  basic_data_phospho_only_helper <- phospho_cln %>%
    anti_join(proteins_cln, 
              by = c("uniprot_acc" = "uniprot_acc",
                     "comparison" = "comparison")) %>%
    anti_join(basic_data_shared,
              by = c("sites_id" = "sites_id",
                     "comparison" = "comparison"))
  
  basic_data_phospho_only <- phospho_tbl %>%
    inner_join(basic_data_phospho_only_helper %>%
                 dplyr::select(phos_row_ids, comparison),
               by = c("phos_row_ids", "comparison")) %>%
    dplyr::mutate(norm_phos_logFC = log2FC, combined_q_mod = q.mod) %>%
    dplyr::rename(log2FC.phos= log2FC, q.mod.phos = q.mod) %>%
    dplyr::mutate(status = "Phos_Only")
  
  if("fdr.mod" %in% colnames(proteins_cln) &
     "fdr.mod" %in% colnames(phospho_cln)) {
    basic_data_phospho_only <- basic_data_phospho_only %>%
      dplyr::rename(fdr.mod.phos = fdr.mod)
  }
  
  list_of_phospho_only_columns <- c("comparison", "norm_phos_logFC", "combined_q_mod", "sites_id", 
                                    "uniprot_acc", "position", "residue", "sequence", "log2FC.phos", 
                                    "q.mod.phos", "fdr.mod.phos", "status", "phos_row_ids" )
  
  if(!"combined_fdr_mod" %in% colnames(basic_data_shared)) {
    list_of_phospho_only_columns <- 
      list_of_phospho_only_columns[!list_of_phospho_only_columns %in% 
                                     c("combined_fdr_mod", "fdr.mod.phos")]
  }
  
  basic_data_phospho_only <- basic_data_phospho_only %>%
    dplyr::select(one_of(list_of_phospho_only_columns)) %>%
    arrange(comparison, combined_q_mod, norm_phos_logFC) %>%
    distinct() %>%
    as.data.frame
  
  basic_data <- basic_data_shared %>%
    bind_rows(basic_data_phospho_only %>%
                dplyr::mutate( position = purrr::map_chr(position, as.character)))
  
  # if (nrow(basic_data %>% distinct(sites_id)) != nrow(phospho_cln %>% distinct(sites_id))){
  #    logwarn("nrow(basic_data) != nrow(phospho_cln)")
  # }
  
  annotation_from_phospho_tbl <- phospho_tbl_orig %>%
    dplyr::mutate(phos_row_ids = sites_id) %>%
    dplyr::select(-q.mod, -p.mod, -log2FC, -uniprot_acc, -position, -residue, -sequence)
  
  annotated_phos_tbl_with_repeats <- basic_data %>%
    left_join(annotation_from_phospho_tbl, 
              by = c("sites_id" = "sites_id", 
                     "comparison" = "comparison",
                     "phos_row_ids" = "phos_row_ids")) %>%
    dplyr::select(!matches("log2norm\\.\\d+\\.(left|right)") & 
                    !matches("raw\\.\\d+\\.(left|right)")) %>%
    arrange(comparison, combined_q_mod, norm_phos_logFC) %>%
    distinct()
  
  vroom::vroom_write(annotated_phos_tbl_with_repeats, 
                     file.path(output_dir, 
                               "norm_phosphosite_lfc_minus_protein_lfc_with_protein_repeats.tsv"))
  
  list_of_long_columns <- intersect(colnames(basic_data),
                                    c("protein_names", "ENSEMBL", "PROTEIN-NAMES", "KEYWORDS",
                                      "GO-ID", "go_biological_process", "go_cellular_compartment",
                                      "go_molecular_function", "reactome_term", "majority_protein_ids"))
  
  writexl::write_xlsx(annotated_phos_tbl_with_repeats %>%
                        mutate(across(any_of(colnames(annotated_phos_tbl_with_repeats)), 
                                      \(x){substr(x, 1, 32760)})),
                      file.path(output_dir,
                                "norm_phosphosite_lfc_minus_protein_lfc_with_protein_repeats.xlsx"))
  
  #----- norm_phosphosite_lfc_minus_protein_lfc_basic_no_repeast
  #- Given an input table with sites_id and uniprot_acc columns, 
  #- work out the ranking of the uniprot_acc within the sites_id
  getUniprotAccRankFromSitesId <- function(input_table, 
                                           uniprot_acc, 
                                           sites_id) {
    input_table %>%
      dplyr::mutate(uniprot_acc_split = str_split({{uniprot_acc}}, ":" ) %>% 
                      purrr::map_chr(1)) %>%
      dplyr::mutate(uniprot_list = str_split({{sites_id}}, "!") %>% 
                      purrr::map_chr(1) %>% 
                      str_split( ":")) %>%
      dplyr::mutate(gene_list_position = purrr::map2_int(uniprot_acc_split, 
                                                         uniprot_list, 
                                                         ~{which(.x == .y)[1]})) %>%
      relocate(uniprot_acc_split, .after=lazyeval::as_name(enquo(sites_id))) %>%
      relocate(uniprot_list, .after="uniprot_acc_split") %>%
      relocate(gene_list_position, .after="uniprot_list")
  }
  
  if(protein_group_filter == "protein_rank") {
    temp_ranking <- getUniprotAccRankFromSitesId(basic_data, uniprot_acc, sites_id)
    min_ranking <- temp_ranking %>%
      group_by(sites_id) %>%
      summarise(gene_list_position = min(gene_list_position)) %>%
      ungroup()
    basic_data_no_repeats <- temp_ranking %>%
      dplyr::inner_join(min_ranking, 
                        by = c("sites_id", "gene_list_position")) %>%
      dplyr::select(-uniprot_list, -gene_list_position, -uniprot_acc_split)
  } else if(protein_group_filter == "best_p_value") {
    basic_data_no_repeats <- getUniprotAccRankFromSitesId(basic_data, uniprot_acc, sites_id) %>%
      inner_join(basic_data %>% 
                   group_by(comparison, sites_id) %>%
                   summarise(combined_q_mod = min(combined_q_mod)) %>%
                   ungroup(),
                 by = c("comparison",
                        "combined_q_mod",
                        "sites_id")) %>%
      group_by(comparison, sites_id) %>%
      mutate(row_rank = dense_rank(gene_list_position)) %>%
      ungroup() %>%
      dplyr::filter(row_rank == 1) %>%
      dplyr::select(-uniprot_list, -gene_list_position, -uniprot_acc_split, -row_rank)
  } else if(protein_group_filter == "best_log_fc") {
    basic_data_no_repeats <- getUniprotAccRankFromSitesId(basic_data, uniprot_acc, sites_id) %>%
      mutate(abs_norm_phos_logFC = abs(norm_phos_logFC)) %>%
      inner_join(basic_data %>%
                   mutate(abs_norm_phos_logFC = abs(norm_phos_logFC)) %>%
                   group_by(comparison, sites_id) %>%
                   summarise(abs_norm_phos_logFC = min(abs_norm_phos_logFC)) %>%
                   ungroup(),
                 by = c("comparison",
                        "sites_id",
                        "abs_norm_phos_logFC")) %>%
      group_by(comparison, sites_id) %>%
      mutate(row_rank = dense_rank(gene_list_position)) %>%
      ungroup() %>%
      dplyr::filter(row_rank == 1) %>%
      dplyr::select(-uniprot_list, -gene_list_position, -uniprot_acc_split, -row_rank, -abs_norm_phos_logFC)
  } else {
    print(paste0("Invalid protein group repeats filter: '", protein_group_filter, "'."))
  }
  
  vroom::vroom_write(basic_data_no_repeats, 
                     file.path(output_dir, 
                               "norm_phosphosite_lfc_minus_protein_lfc_basic_no_repeast.tsv"))
  
  #----- Join normalized phosphopeptide abundance table with phosphosite annotations
  annotation_from_phospho_tbl <- phospho_tbl_orig %>%
    dplyr::mutate(phos_row_ids = sites_id) %>%
    dplyr::select(-q.mod, -p.mod, -log2FC, -uniprot_acc,  -position, -residue, -sequence)
  
  annotated_phos_tbl <- basic_data_no_repeats %>%
    left_join(annotation_from_phospho_tbl, 
              by = c("sites_id" = "sites_id",
                     "comparison" = "comparison",
                     "phos_row_ids" = "phos_row_ids")) %>%
    dplyr::select(!matches( "log2norm\\.\\d+\\.(left|right)") & 
                    !matches("raw\\.\\d+\\.(left|right)")) %>%
    arrange(comparison, combined_q_mod, norm_phos_logFC) %>%
    distinct()
  
  vroom::vroom_write(annotated_phos_tbl,
                     file.path(output_dir,
                               "norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv"))
  
  list_of_long_columns <- intersect(colnames(annotated_phos_tbl),
                                    c("protein_names", "ENSEMBL", "PROTEIN-NAMES", "KEYWORDS",
                                      "GO-ID", "go_biological_process", "go_cellular_compartment",
                                      "go_molecular_function", "reactome_term", "majority_protein_ids"))
  
  writexl::write_xlsx(annotated_phos_tbl %>%
                        mutate(across(any_of(colnames(annotated_phos_tbl)), \(x){substr(x, 1, 32760)})),
                      file.path(output_dir,
                                "norm_phosphosite_lfc_minus_protein_lfc_annotated.xlsx"))
  
  #----- Compare before and after normalization with protein abundance
  before_prot_norm <- phospho_cln %>%
    dplyr::filter(q.mod < 0.05) %>%
    dplyr::distinct(comparison, sites_id) %>%
    dplyr::mutate(Normalization = "Before")
  
  after_prot_norm <- basic_data_no_repeats %>%
    dplyr::filter(combined_q_mod < 0.05) %>%
    dplyr::distinct(comparison, sites_id) %>%
    dplyr::mutate(Normalization = "After")
  
  comparisons_order <- before_prot_norm %>%
    distinct(comparison) %>%
    pull(comparison)
  
  compare_before_and_after <- before_prot_norm %>%
    full_join(after_prot_norm, 
              by = c("comparison", "sites_id"), 
              suffix = c(".before", ".after")) %>%
    mutate(Normalization = case_when(is.na(Normalization.before) & !is.na(Normalization.after) ~ "After Only",
                                     !is.na(Normalization.before) & !is.na(Normalization.after) ~ "Before & After",
                                     !is.na(Normalization.before) & is.na(Normalization.after) ~ "Before Only",
                                     TRUE ~ NA_character_ )) %>%
    mutate(Normalization = factor(Normalization,
                                  levels = c("Before Only",
                                             "Before & After",
                                             "After Only"))) %>%
    mutate(comparison = factor(comparison, levels = comparisons_order)) %>%
    group_by(comparison, Normalization) %>%
    summarise(Counts = n()) %>%
    ungroup()
  
  vroom::vroom_write(compare_before_and_after,
                     file.path(output_dir, 
                               "compare_before_and_after_norm_by_prot_abundance.tsv"))
  
  cmp_before_after_plot <- compare_before_and_after %>%
    ggplot(aes(Normalization, Counts)) +
    geom_col() +
    facet_grid(. ~ comparison) +
    geom_text(stat='identity', aes(label= Counts), vjust=-0.5) +
    theme(axis.text.x = element_text(angle = 90))
  
  for(format_ext in plots_format) {
    file_name<-file.path(output_dir,
                         paste0("compare_before_and_after_norm_by_prot_abundance.", format_ext))
    ggsave(plot = cmp_before_after_plot, file_name, width = 14, height = 10)
  }
  
  #----- Count num. differentiall abundant phosphosites after normalization by protein abundance
  #- Count the number of statistically significant differentially expressed proteins 
  #- (according to user-defined threshold)
  countStatDeGenes <- function(data,
                               lfc_thresh = 0,
                               q_val_thresh = 0.05,
                               log_fc_column = log2FC,
                               q_value_column = q.mod) {
    # comparison <- as.data.frame(data) |>
    # distinct(comparison) |>
    # pull(comparison)
    selected_data <- data |>
      dplyr::mutate(status = case_when({{q_value_column}} >= q_val_thresh ~ "Not significant",
                                       {{log_fc_column}} >= lfc_thresh & 
                                         {{q_value_column}} < q_val_thresh ~ "Significant and Up",
                                       {{log_fc_column}} < lfc_thresh & 
                                         {{q_value_column}} < q_val_thresh ~ "Significant and Down",
                                       TRUE ~ "Not significant"))
    counts <- selected_data |>
      group_by(status) |>
      summarise(counts = n()) |>
      ungroup()
    all_possible_status <- data.frame(status = c("Not significant", 
                                                 "Significant and Up", 
                                                 "Significant and Down"))
    results <- all_possible_status |>
      left_join(counts, by = c("status" = "status")) |>
      mutate(counts = ifelse(is.na(counts), 0, counts))
    return(results)
  }
  
  output_counts_tbl <- annotated_phos_tbl %>%
    dplyr::select(comparison, norm_phos_logFC, combined_q_mod) %>%
    group_by(comparison) %>%
    nest(data = c(norm_phos_logFC, combined_q_mod)) %>%
    ungroup %>%
    mutate(counts = purrr::map(data, function(x) { 
      countStatDeGenes(x, 
                       lfc_thresh = 0,
                       q_val_thresh = 0.05,
                       log_fc_column = norm_phos_logFC,
                       q_value_column = combined_q_mod )
    })) %>%
    dplyr::select(-data) %>%
    unnest(counts)
  
  vroom::vroom_write(output_counts_tbl,
                     file.path(output_dir,
                               "num_sig_diff_abundant_norm_phosphosites.tab"))
  
  #----- Create Volcano Plot
  qm.threshold <- 0.05
  logFC.threshold <- 1
  
  basic_data_volcano_plot_data <- basic_data_no_repeats %>%
    left_join(annotated_phos_tbl %>%
                dplyr::distinct(uniprot_acc, gene_name),
              by = c("uniprot_acc" = "uniprot_acc")) %>%
    dplyr::mutate(colour = case_when(abs(norm_phos_logFC) >= logFC.threshold & 
                                       combined_q_mod >= qm.threshold ~ "orange",
                                     abs(norm_phos_logFC) >= logFC.threshold & 
                                       combined_q_mod < qm.threshold ~ "purple",
                                     abs(norm_phos_logFC) < logFC.threshold & 
                                       combined_q_mod < qm.threshold ~ "blue",
                                     TRUE ~ "black")) %>%
    dplyr::mutate(colour = factor(colour, levels=c("orange", "purple", "blue" ,"black")))
  
  basic_data_volcano_plot <- basic_data_volcano_plot_data %>%
    ggplot(aes(norm_phos_logFC, -log10(combined_q_mod), col=colour, key=gene_name)) +
    geom_point() +
    facet_grid(. ~ comparison) +
    scale_colour_manual(values = c("orange", "purple", "blue" ,"black"),
                        labels = c(paste0("Not significant, logFC > ", logFC.threshold),
                                   paste0("Significant, logFC >= ", logFC.threshold),
                                   paste0("Significant, logFC <", logFC.threshold),
                                   "Not Significant"))
  
  for(format_ext in plots_format) {
    file_name <- file.path(output_dir, 
                           paste0("volplot_gg_phos_vs_prot_all.", format_ext))
    ggsave(filename = file_name, plot = basic_data_volcano_plot, width = 15, height = 6)
  }
  
  # ggplotly(basic_data_volcano_plot)
}