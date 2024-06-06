#---------- Phosphoproteome Pathways
#----- Packages
# suppressPackageStartupMessages({
#   p_load(clusterProfiler)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for temporary files
# tmp_dir <- "results/proteomics/cache"
# #- Results table in which the phosphorylation log fold-change is normalized by protein log fold-change
# norm_phos_logfc_file <- "results/phosphoproteomics/norm_phos_by_prot_abundance/norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv"
# #- Column name in the input file that contains the log fold-change values of the phosphosites
# log_fc_column_name <- "norm_phos_logFC"
# #- Column name in the input file that contains the false discovery rate values of the phosphosites
# fdr_column_name <- "combined_q_mod"
# #- Column name in the input file that contains the gene names
# results_gene_name_column_name <- "gene_name"
# #- File with table of differentially expressed proteins
# proteins_file <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
# #- phosphoproteins or proteins_and_phosphoproteins
# background <- "proteins_and_phosphoproteins"
# #- The protein id in the annotation file
# protein_id <- "uniprot_acc"
# #- File that converts annotation ID to annotation name
# dictionary_file <- "data/HomoSapiens/uniprot/go_terms_table_python_all.tab"
# #- The maximum number of genes associaed with each gene set
# max_gene_set_size <- list("250,500")
# #- The minimum number of genes associaed with each gene set
# min_gene_set_size <- list("5,10")
# #- p-value threshold below which a phosphosite is significantly enriched
# site_p_val_thresh <- 0.05
# #- p-value threshold below which a GO term is significantly enriched
# p_val_thresh <- 0.05
# #- A file that contains the dictionary to convert uniprot accession to gene symbol
# uniprot_to_gene_symbol_file <- "data/HomoSapiens/uniprot/data.tab"
# #- The name of the column that contained the protein ID to convert into gene symbol
# protein_id_lookup_column <- "Entry"
# #- The name of the column that contained the gene symbol
# gene_symbol_column <- "Gene Names"

#----- GO
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phosphoproteins_go"
# #- File with protein accession and functional annotation data
# annotation_file <- "data/HomoSapiens/uniprot/go_terms_table_python_all.tab"
# #- Name of the functional annotation data, added to the output file name
# annotation_type <- "go_enrichment"
# #- The annotatio id in the annotation file
# annotation_id <- "go_id"
# #- Column with the long name and biological details of the functional annotation
# annotation_column <- "go_term" 
# #- The aspect of the GO term
# aspect_column <- "go_type"

#----- Reactome
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phosphoproteins_reactome"
# #- File with protein accession and functional annotation data
# annotation_file <- "data/HomoSapiens/reactome/reactome_data_with_header.txt"
# #- File that converts annotation ID to annotation name
# dictionary_file <- "data/HomoSapiens/reactome/reactome_data_with_header.txt"
# #- Name of the functional annotation data, added to the output file name
# annotation_type <- "reactome_enrichment"
# #- The annotatio id in the annotation file
# annotation_id <- "reactome_id"
# #- Column with the long name and biological details of the functional annotation
# annotation_column <- "reactome_term"
# #- The aspect of the GO term
# aspect_column <- NULL

#----- Kegg
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phosphoproteins_kegg"
# #- File with protein accession and functional annotation data
# annotation_file <- "data/HomoSapiens/kegg/gene_sets.tab"
# #- File that converts annotation ID to annotation name
# dictionary_file <- "data/HomoSapiens/kegg/gene_sets.tab"
# #- Name of the functional annotation data, added to the output file name
# annotation_type <- "kegg_enrichment"
# #- The annotatio id in the annotation file
# annotation_id <- "pathway_id"
# #- Column with the long name and biological details of the functional annotation
# annotation_column <- "pathway_name"
# #- The aspect of the GO term
# aspect_column <- NULL

#----- Function
getPhosphoproteomePathways <- function(
    tmp_dir,
    output_dir,
    norm_phos_logfc_file,
    log_fc_column_name,
    fdr_column_name,
    results_gene_name_column_name,
    proteins_file,
    background,
    annotation_file,
    protein_id,
    annotation_id,
    annotation_column,
    annotation_type,
    aspect_column,
    dictionary_file,
    max_gene_set_size,
    min_gene_set_size,
    site_p_val_thresh,
    p_val_thresh,
    uniprot_to_gene_symbol_file,
    protein_id_lookup_column,
    gene_symbol_column
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
  
  #----- Read phosphosite log fold-change normalized by protein log fold-change file
  de_phos <- vroom::vroom(norm_phos_logfc_file)
  
  #----- Check log fold-change and FDR column names exists in input file
  # if (!log_fc_column_name %in% colnames(de_phos)) {
  #   print("Column '%s' is not found in the input table.", log_fc_column_name )
  # }
  # if (!fdr_column_name %in% colnames(de_phos)) {
  #   print("Column '%s' is not found in the input table.", fdr_column_name)
  # }
  # if(!log_fc_column_name %in% colnames(de_phos) | 
  #    !fdr_column_name %in% colnames(de_phos)) {
  #   stop()
  # }
  
  #----- Positive phosphoproteins
  positive_phosphoproteins <- de_phos %>%
    dplyr::filter(!!rlang::sym(fdr_column_name) < site_p_val_thresh & 
                    !!rlang::sym(log_fc_column_name) > 0) %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    #- Strip away isoform information
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
    mutate(gene_name_first = purrr::map_chr(!!rlang::sym(results_gene_name_column_name), 
                                            ~str_split(., ":") %>% 
                                              map_chr(1))) %>%
    mutate(protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, 
                                               ~str_split(., ":") %>% 
                                                 map_chr(1))) %>%
    group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
    summarise(max_norm_phos_logFC = max(!!rlang::sym(log_fc_column_name))) %>%
    ungroup() %>%
    arrange(comparison, desc(max_norm_phos_logFC))
  
  vroom::vroom_write(positive_phosphoproteins,
                     file.path(output_dir,
                               "all_phosphoproteins_with_positive_logFC_sites.tab"))
  
  list_of_comparisons <- positive_phosphoproteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) { 
    positive_phosphoproteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "all_phosphoproteins_with_positive_logFC_sites.tab"),
                         col_names=FALSE)
  })
  
  #----- Negative phosphoproteins
  negative_phosphoproteins <- de_phos %>%
    dplyr::filter(!!rlang::sym(fdr_column_name) < site_p_val_thresh & 
                    !!rlang::sym(log_fc_column_name) < 0) %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    #- Strip away isoform information
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
    mutate(gene_name_first = purrr::map_chr(!!rlang::sym(results_gene_name_column_name), 
                                            ~str_split(., ":") %>% 
                                              map_chr(1))) %>%
    mutate(protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, 
                                               ~str_split(., ":") %>% 
                                                 map_chr(1))) %>%
    group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
    summarise(min_norm_phos_logFC = min(!!rlang::sym(log_fc_column_name))) %>%
    ungroup() %>%
    arrange(comparison, min_norm_phos_logFC)
  
  vroom::vroom_write(negative_phosphoproteins,
                     file.path(output_dir,
                               "all_phosphoproteins_with_negative_logFC_sites.tab"))
  
  list_of_comparisons <- negative_phosphoproteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    negative_phosphoproteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "all_phosphoproteins_with_negative_logFC_sites.tab"),
                         col_names=FALSE)
  })
  
  #----- Positive Only Phosphoproteins
  positive_only_phosphoproteins <- positive_phosphoproteins %>%
    anti_join(negative_phosphoproteins, 
              by = c("uniprot_acc_first" = "uniprot_acc_first"))
  
  vroom::vroom_write(positive_only_phosphoproteins,
                     file.path(output_dir,
                               "phosphoproteins_with_only_positive_logFC_sites.tab"))
  
  list_of_comparisons <- positive_only_phosphoproteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    positive_only_phosphoproteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_only_positive_logFC_sites.tab"),
                         col_names=FALSE)
  })
  
  #----- Negative Only Phosphoproteins
  negative_only_phosphoproteins <- negative_phosphoproteins %>%
    anti_join(positive_phosphoproteins, 
              by = c("uniprot_acc_first" = "uniprot_acc_first"))
  
  vroom::vroom_write(negative_only_phosphoproteins,
                     file.path(output_dir,
                               "phosphoproteins_with_only_negative_logFC_sites.tab"))
  
  list_of_comparisons <- negative_only_phosphoproteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    negative_only_phosphoproteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_only_negative_logFC_sites.tab"),
                         col_names=FALSE)
  })
  
  #----- Overlapping Phosphoproteins
  overlapping_phosphoproteins <- positive_phosphoproteins %>%
    inner_join(negative_phosphoproteins, 
               by=c("uniprot_acc_first" = "uniprot_acc_first",
                    "comparison" = "comparison",
                    "gene_name_first" = "gene_name_first",
                    "protein_name_first" = "protein_name_first"))
  
  vroom::vroom_write(overlapping_phosphoproteins,
                     file.path(output_dir,
                               "phosphoproteins_with_positive_and_negative_logFC_sites.tab"))
  
  list_of_comparisons <- overlapping_phosphoproteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison){
    overlapping_phosphoproteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "phosphoproteins_with_positive_and_negative_logFC_sites.tab"),
                         col_names=FALSE)
  })
  
  #----- All Phosphoproteins with significant sites
  all_phosphoproteins_with_significant_da_sites <- positive_phosphoproteins %>%
    bind_rows(negative_phosphoproteins) %>%
    distinct()
  
  vroom::vroom_write(all_phosphoproteins_with_significant_da_sites,
                     file.path(output_dir,
                               "all_phosphoproteins_with_significant_da_sites.tab"))
  
  list_of_comparisons <- all_phosphoproteins_with_significant_da_sites %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    all_phosphoproteins_with_significant_da_sites %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "all_phosphoproteins_with_significant_da_sites.tab"),
                         col_names=FALSE)
  })
  
  #- Classify phosphoproteins as up and down by total log2 FC of significant phosphosites
  group_phosphoproteins_by_phosphosites_lfc_total <- de_phos  %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>%
    dplyr::filter(!!rlang::sym(fdr_column_name) < site_p_val_thresh) %>%
    group_by(comparison, uniprot_acc_first) %>%
    summarise(total_log2FC = sum(!!rlang::sym(log_fc_column_name))) %>%
    ungroup() %>%
    mutate(direction = case_when(total_log2FC > 0 ~ "up",
                                 total_log2FC < 0 ~ "down",
                                 TRUE ~ "no_change")) %>%
    dplyr::filter(direction != "no_change")
  
  gene_names_list <- de_phos %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    #- Strip away isoform information
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
    mutate(gene_name_first = purrr::map_chr(!!rlang::sym(results_gene_name_column_name), 
                                            ~str_split(., ":") %>% 
                                              map_chr(1))) %>%
    mutate( protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, 
                                                ~str_split(., ":") %>% 
                                                  map_chr(1))) %>%
    distinct(comparison, uniprot_acc_first, gene_name_first, protein_name_first)
  
  group_phosphoproteins_by_phosphosites_lfc_total_up <- group_phosphoproteins_by_phosphosites_lfc_total %>%
    dplyr::filter(direction == "up") %>%
    left_join(gene_names_list, 
              by = c("uniprot_acc_first" = "uniprot_acc_first",
                     "comparison" = "comparison"))
  
  vroom::vroom_write(group_phosphoproteins_by_phosphosites_lfc_total_up,
                     file.path(output_dir,
                               "group_phosphoproteins_by_phosphosites_lfc_total_up.tab"))
  
  list_of_comparisons <- group_phosphoproteins_by_phosphosites_lfc_total_up %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    group_phosphoproteins_by_phosphosites_lfc_total_up %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "group_phosphoproteins_by_phosphosites_lfc_total_up.tab"),
                         col_names=FALSE)
  })
  
  group_phosphoproteins_by_phosphosites_lfc_total_down <- group_phosphoproteins_by_phosphosites_lfc_total %>%
    dplyr::filter(direction == "down") %>%
    left_join(gene_names_list, 
              by = c("uniprot_acc_first" = "uniprot_acc_first",
                     "comparison" = "comparison"))
  
  vroom::vroom_write(group_phosphoproteins_by_phosphosites_lfc_total_down,
                     file.path(output_dir,
                               "group_phosphoproteins_by_phosphosites_lfc_total_down.tab"))
  
  list_of_comparisons <- group_phosphoproteins_by_phosphosites_lfc_total_down %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, 
              ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    group_phosphoproteins_by_phosphosites_lfc_total_down %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "group_phosphoproteins_by_phosphosites_lfc_total_down.tab" ),
                         col_names=FALSE)
  })
  
  #----- Read proteomics abundance data
  proteins_tbl_orig <- vroom::vroom(proteins_file, delim="\t")
  
  #----- Background Phosphoproteins
  background_phosphoproteins  <- de_phos %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    #- Strip away isoform information
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
    mutate(gene_name_first = purrr::map_chr(!!rlang::sym(results_gene_name_column_name), 
                                            ~str_split(., ":") %>% 
                                              map_chr(1))) %>%
    mutate(protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, 
                                               ~str_split(., ":") %>% 
                                                 map_chr(1))) %>%
    distinct(uniprot_acc_first) %>%
    arrange(uniprot_acc_first)
  
  vroom::vroom_write(background_phosphoproteins,
                     file.path(output_dir, "background_phosphoproteins.tab"),
                     col_names=FALSE)
  
  #----- Background Proteins
  background_proteins <- proteins_tbl_orig %>%
    distinct(uniprot_acc) %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, 
                                              ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>%
    distinct(uniprot_acc_first)
  
  vroom::vroom_write(background_proteins,
                     file.path(output_dir, "background_proteins.tab"),
                     col_names=FALSE)
  
  #----- Background Proteins Phosphoproteins
  background_proteins_phosphoproteins <- background_proteins %>%
    bind_rows(background_phosphoproteins) %>%
    distinct(uniprot_acc_first)
  
  vroom::vroom_write(background_proteins_phosphoproteins,
                     file.path(output_dir, "background_proteins_phosphoproteins.tab" ),
                     col_names=FALSE )
  
  #----- Pathways
  if(is.null(annotation_file)) {
    stop("No annotation file provided.")
  }
  
  buildAnnotationIdToAnnotationNameDictionary <- function(input_table, 
                                                          annotation_column, 
                                                          annotation_id_column) {
    id_to_annotation_dictionary <- NA
    dictionary_pair <- input_table %>%
      dplyr::filter( !is.na({{annotation_column}}) & !is.na({{annotation_id_column}})) %>%
      distinct({{annotation_column}},
               {{annotation_id_column}})
    id_to_annotation_dictionary <- dictionary_pair %>%
      pull({{annotation_column}} )
    names(id_to_annotation_dictionary ) <- dictionary_pair %>%
      pull( {{annotation_id_column}})
    id_to_annotation_dictionary
  }
  
  #- Tidy up the annotation ID to annotation term name dictionary
  dictionary <- vroom::vroom(dictionary_file)
  
  id_to_annotation_dictionary <- 
    buildAnnotationIdToAnnotationNameDictionary(input_table = dictionary,
                                                annotation_column = !!rlang::sym(annotation_column),
                                                annotation_id_column = !!rlang::sym(annotation_id))
  
  #- preparing the enrichment test
  go_annot <- vroom::vroom(annotation_file)
  
  background_list <- background_phosphoproteins
  if(background == "phosphoproteins") {
    background_list <- background_phosphoproteins
  } else if(background == "proteins_and_phosphoproteins") {
    background_list <- background_proteins_phosphoproteins
  }
  
  #- Function used for parsing a list of minimum or maximum gene set size from command line
  parseNumList <- function(input_text) {
    if(str_detect(input_text, "[.,;:]")) {
      str_split(input_text, "[.,;:]")[[1]] %>%
        purrr::map_int(as.integer)
    } else {
      return(as.integer(input_text))
    }
  }
  
  min_gene_set_size_list <- parseNumList(min_gene_set_size)
  max_gene_set_size_list <- parseNumList(max_gene_set_size)
  
  list_of_comparisons <- all_phosphoproteins_with_significant_da_sites %>%
    distinct(comparison) %>%
    pull(comparison)
  
  #- Tidy up GO aspect list, marked as null if not using GO terms
  go_aspect_list <- NA
  if(!is.null(aspect_column)) {
    go_aspect_list <- go_annot %>%
      dplyr::filter(!is.na(!!rlang::sym(aspect_column))) %>%
      distinct(!!rlang::sym(aspect_column)) %>%
      pull(!!rlang::sym(aspect_column) ) # c("C", "F", "P")
  } else {
    go_aspect_list <- NA
  }
  
  list_of_genes_list <- list(all_significant = all_phosphoproteins_with_significant_da_sites,
                             overlap_only = overlapping_phosphoproteins,
                             negative_only = negative_only_phosphoproteins,
                             positive_only = positive_only_phosphoproteins,
                             negative_plus_overlap = negative_phosphoproteins,
                             positive_plus_overlap = positive_phosphoproteins,
                             negative_sum_sig_phosphosites = group_phosphoproteins_by_phosphosites_lfc_total_down,
                             positive_sum_sig_phosphosites = group_phosphoproteins_by_phosphosites_lfc_total_up)
  
  input_params <- expand_grid(names_of_genes_list = names(list_of_genes_list),
                              go_aspect = go_aspect_list,
                              input_comparison = list_of_comparisons,
                              min_size = min_gene_set_size_list,
                              max_size = max_gene_set_size_list)
  
  input_params_updated <- input_params %>%
    mutate(input_table = purrr::map(names_of_genes_list, \(x) list_of_genes_list[[x]]))
  
  #----- Run Enrichment
  convertIdToAnnotation <- function(id, id_to_annotation_dictionary) {
    return(ifelse(!is.null(id_to_annotation_dictionary[[id]]), 
                  id_to_annotation_dictionary[[id]], 
                  NA_character_))
  }
  
  oneGoEnrichment <- function(go_annot, 
                              background_list, 
                              go_aspect, 
                              query_list, 
                              id_to_annotation_dictionary,
                              annotation_id, 
                              protein_id, 
                              aspect_column, 
                              p_val_thresh, 
                              min_gene_set_size, 
                              max_gene_set_size) {
    join_condition <- rlang::set_names(c(colnames(background_list)[1]),
                                       c(as_name(enquo(protein_id))))
    if(!is.na(go_aspect)) {
      go_annot_filt <- go_annot %>%
        dplyr::filter({{aspect_column}} == go_aspect) |>
        mutate({{protein_id}} := purrr::map_chr({{protein_id}}, as.character))
    } else {
      go_annot_filt <- go_annot |>
        mutate({{protein_id}} := purrr::map_chr({{protein_id}}, as.character))
    }
    
    filtered_go_terms <- go_annot_filt %>%
      inner_join(background_list, by = join_condition) %>%
      group_by({{annotation_id}}) %>%
      summarise(counts =n()) %>%
      ungroup() %>%
      arrange(desc(counts)) %>%
      dplyr::filter(counts <= max_gene_set_size & counts >= min_gene_set_size) %>%
      dplyr::select(-counts)
    
    term_to_gene_tbl_filt <- go_annot_filt %>%
      inner_join(background_list, by = join_condition) %>%
      dplyr::inner_join(filtered_go_terms, by = as_name(enquo(annotation_id))) %>%
      dplyr::rename(gene = as_name(enquo(protein_id)),
                    term = as_name(enquo(annotation_id))) %>%
      dplyr::select(term, gene) %>%
      dplyr::distinct(term, gene)
    
    #- Avoid singleton GO terms
    terms_to_avoid <- term_to_gene_tbl_filt %>%
      distinct() %>%
      dplyr::inner_join(data.frame(uniprot_acc = query_list), by=c("gene" = "uniprot_acc")) %>%
      distinct() %>%
      group_by(term) %>%
      summarise(counts =n()) %>%
      ungroup() %>%
      dplyr::filter(counts < 2)
    
    term_to_gene_tbl_filt_no_singleton <- term_to_gene_tbl_filt %>%
      dplyr::anti_join(terms_to_avoid, by="term")
    
    no_singleton_terms_query_gene_list <- intersect(query_list,
                                                    term_to_gene_tbl_filt_no_singleton %>%
                                                      dplyr::distinct(gene) %>%
                                                      dplyr::pull(gene))
    
    enrichment_result <- enricher(no_singleton_terms_query_gene_list,
                                  pvalueCutoff = p_val_thresh, 
                                  pAdjustMethod = "BH", 
                                  minGSSize = min_gene_set_size, 
                                  maxGSSize = max_gene_set_size, 
                                  qvalueCutoff = p_val_thresh, 
                                  TERM2GENE =term_to_gene_tbl_filt_no_singleton)
    
    if(!is.null(enrichment_result)) {
      output_table <- as.data.frame(enrichment_result) %>%
        dplyr::mutate(term = purrr::map_chr(ID, function(id) {
          convertIdToAnnotation(id, id_to_annotation_dictionary)})) %>%
        dplyr::relocate(term, .before="Description") %>%
        dplyr::mutate(min_gene_set_size = min_gene_set_size,
                      max_gene_set_size = max_gene_set_size )
      
      output_table_with_go_aspect <- NA
      if(!is.na(go_aspect)) {
        output_table_with_go_aspect <- output_table %>%
          dplyr::mutate({{aspect_column}} := go_aspect)
      } else {
        output_table_with_go_aspect <- output_table
      }
      return(output_table_with_go_aspect)
    } else {
      return( NULL)
    }
  }
  
  runOneGoEnrichmentInOutFunction <- function(comparison_column,
                                              protein_id_column,
                                              go_annot = go_annot,
                                              background_list,
                                              id_to_annotation_dictionary,
                                              annotation_id,
                                              protein_id,
                                              aspect_column,
                                              p_val_thresh,
                                              names_of_genes_list,
                                              input_table,
                                              go_aspect,
                                              input_comparison,
                                              min_gene_set_size,
                                              max_gene_set_size) {
    oneGoEnrichmentPartial <- NA
    if(!is.null(aspect_column)) {
      oneGoEnrichmentPartial <- purrr::partial(oneGoEnrichment,
                                               go_annot = go_annot,
                                               background_list = background_list,
                                               id_to_annotation_dictionary = id_to_annotation_dictionary,
                                               annotation_id = {{annotation_id}},
                                               protein_id = {{protein_id}},
                                               aspect_column = !!rlang::sym(aspect_column),
                                               p_val_thresh = p_val_thresh)
    } else {
      oneGoEnrichmentPartial <- purrr::partial(oneGoEnrichment,
                                               go_annot = go_annot,
                                               background_list = background_list,
                                               id_to_annotation_dictionary = id_to_annotation_dictionary,
                                               annotation_id = {{annotation_id}},
                                               protein_id = {{protein_id}},
                                               aspect_column = aspect_column,
                                               p_val_thresh = p_val_thresh)
    }
    query_list <- input_table %>%
      dplyr::filter({{comparison_column}} == input_comparison) %>%
      pull({{protein_id_column}})
    
    enrichment_temp <- oneGoEnrichmentPartial(go_aspect = go_aspect, 
                                              query_list = query_list, 
                                              min_gene_set_size = min_gene_set_size, 
                                              max_gene_set_size = max_gene_set_size)
    if(!is.null(enrichment_temp)) {
      enrichment_result <- enrichment_temp %>%
        dplyr::mutate({{comparison_column}} := input_comparison) %>%
        dplyr::mutate(names_of_genes_list = names_of_genes_list)
      return(enrichment_result )
    } else {
      return (NULL)
    }
  }
  
  runOneGoEnrichmentInOutFunctionPartial <- purrr::partial(runOneGoEnrichmentInOutFunction,
                                                           comparison_column = comparison,
                                                           protein_id_column = uniprot_acc_first,
                                                           go_annot = go_annot,
                                                           background_list = background_list,
                                                           id_to_annotation_dictionary = id_to_annotation_dictionary,
                                                           annotation_id = !!rlang::sym(annotation_id),
                                                           protein_id = !!rlang::sym(protein_id),
                                                           aspect_column = aspect_column,
                                                           p_val_thresh = p_val_thresh)
  
  enrichment_result <- input_params_updated %>%
    mutate(enrichment_results = purrr::pmap(list(names_of_genes_list, 
                                                 input_table, 
                                                 go_aspect, 
                                                 input_comparison, 
                                                 min_size, 
                                                 max_size) ,
                                            \(names_of_genes_list, 
                                              input_table, 
                                              go_aspect, 
                                              input_comparison, 
                                              min_size, 
                                              max_size) {runOneGoEnrichmentInOutFunctionPartial(
                                                names_of_genes_list = names_of_genes_list,
                                                input_table = input_table,
                                                go_aspect = go_aspect,
                                                input_comparison = input_comparison,
                                                min_gene_set_size = min_size,
                                                max_gene_set_size = max_size)})) %>%
    dplyr::select(enrichment_results) %>%
    unnest()
  
  getUniprotAccToGeneSymbolDictionary <- function(input_table,
                                                  protein_id_lookup_column,
                                                  gene_symbol_column,
                                                  protein_id) {
    #- Clean up protein ID to gene sybmol table
    uniprot_to_gene_symbol <- input_table %>%
      dplyr::select({{protein_id_lookup_column}},
                    {{gene_symbol_column}}) %>%
      dplyr::rename({{protein_id}} := as_name(enquo(protein_id_lookup_column))) %>%
      dplyr::rename(gene_symbol = as_name(enquo(gene_symbol_column))) %>%
      dplyr::mutate(gene_symbol = str_split(gene_symbol, " ") %>%
                      purrr::map_chr(1)) %>%
      dplyr::distinct({{protein_id}}, gene_symbol)
    #- Convert to lookup dictionary
    uniprot_to_gene_symbol_dict <- uniprot_to_gene_symbol %>%
      pull(gene_symbol)
    names(uniprot_to_gene_symbol_dict) <- uniprot_to_gene_symbol %>%
      pull({{protein_id}})
    uniprot_to_gene_symbol_dict
  }
  
  convertProteinAccToGeneSymbol <- function(gene_id_list, dictionary) { 
    purrr::map_chr(gene_id_list,
                   ~{ifelse( . %in% names(dictionary ), 
                             dictionary[[.]], 
                             NA_character_)}) %>% 
      paste( collapse="/")
  }
  
  if(is.null(enrichment_result) |
     nrow(enrichment_result) == 0) {
    warnings("No enriched terms were identified.")
  } else {
    enrichment_result_add_gene_symbol <- NA
    #- Convert Uniprot accession to gene names
    if(file.exists(uniprot_to_gene_symbol_file)) {
      # args$uniprot_to_gene_symbol_file <- "/home/ubuntu/Workings/2021/ALPK1_BMP_06/Data/UniProt/data.tab"
      # args$protein_id_lookup_column <- "Entry"
      # args$gene_symbol_column <- "Gene names"
      #- Clean up protein ID to gene symbol table
      uniprot_tab_delimited_tbl <- vroom::vroom(file.path(uniprot_to_gene_symbol_file))
      uniprot_to_gene_symbol_dict <- getUniprotAccToGeneSymbolDictionary(uniprot_tab_delimited_tbl,
                                                                         !!rlang::sym(protein_id_lookup_column),
                                                                         !!rlang::sym(gene_symbol_column),
                                                                         !!rlang::sym(protein_id))
      convertProteinAccToGeneSymbolPartial <- purrr::partial(convertProteinAccToGeneSymbol,
                                                             dictionary = uniprot_to_gene_symbol_dict)
      print(head(enrichment_result))
      enrichment_result_add_gene_symbol <- enrichment_result %>%
        mutate(gene_id_list = str_split(geneID, "/")) %>%
        mutate(gene_symbol = purrr::map_chr(gene_id_list,
                                            convertProteinAccToGeneSymbolPartial)) %>%
        dplyr::select(-gene_id_list)
    } else {
      enrichment_result_add_gene_symbol <- enrichment_result
    }
    ## generate the output files
    purrr::walk(list_of_comparisons, function(input_comparison) {
      output_file <- paste0(annotation_type, "_table_", input_comparison, ".tab")
      vroom::vroom_write(enrichment_result_add_gene_symbol %>%
                           dplyr::filter(comparison == input_comparison),
                         file = file.path(output_dir,
                                          input_comparison,
                                          output_file))})
  }
}