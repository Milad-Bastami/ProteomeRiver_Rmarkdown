#---------- Proteome Pathways
#----- Packages
# suppressPackageStartupMessages({
#   p_load(clusterProfiler)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for temporary files
# tmp_dir <- "results/proteomics/cache"
# #- File with table of differentially expressed proteins
# proteins_file <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
# #- Column name in the input file that contains the log fold-change values of the proteins
# log_fc_column_name <- "log2FC"
# #- Column name in the input file that contains the false discovery rate values of the proteins
# fdr_column_name <- "q.mod"
# #- p-value threshold below which a protein is significantly enriched
# protein_p_val_thresh <- 0.05
# #- The protein id in the annotation file
# protein_id <- "uniprot_acc"
# #- The maximum number of genes associaed with each gene set
# max_gene_set_size <- list("250,500")
# #- The minimum number of genes associaed with each gene set
# min_gene_set_size <- list("5,10")
# #- p-value threshold below which a GO term is significantly enriched
# p_val_thresh <- 0.05
# #- A file that contains the dictionary to convert uniprot accession to gene symbol. 
# #- Uses the column specified in 'protein_id_lookup_column' flag for protein ID. 
# #- Uses the column specified in 'gene_symbol_column' for the gene symbol column
# uniprot_to_gene_symbol_file <- "data/HomoSapiens/uniprot/data.tab"
# #- The name of the column that contained the protein ID to convert into gene symobl
# protein_id_lookup_column <- "Entry"
# #- The name of the column that contained the gene symobl
# gene_symbol_column <- "Gene Names"
# #- A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg]
# plots_format <- list("pdf","png")
# #- Directory path for all results files
# output_dir <- "results/proteomics/de_proteins_go_list"
# #- File with protein accession and functional annotation data
# annotation_file <- "data/HomoSapiens/uniprot/go_terms_table_python_all.tab"
# #- File that converts annotation ID to annotation name
# dictionary_file <- "data/HomoSapiens/uniprot/go_terms_table_python_all.tab"
# #- Name of the functional annotation data, added to the output file name
# annotation_type <- "go_enrichment"
# #- The annotatio id in the annotation file
# annotation_id <- "go_id"
# #- Column with the long name and biological details of the functional annotation
# annotation_column <- "go_term"
# #- The aspect of the GO term
# aspect_column <- "go_type"

#----- Function
getProteomePathways <- function(
    tmp_dir,
    proteins_file,
    log_fc_column_name,
    fdr_column_name,
    protein_p_val_thresh,
    protein_id,
    p_val_thresh,
    uniprot_to_gene_symbol_file,
    protein_id_lookup_column,
    gene_symbol_column,
    max_gene_set_size,
    min_gene_set_size,
    plots_format,
    output_dir,
    annotation_file,
    dictionary_file,
    annotation_type,
    annotation_id,
    annotation_column,
    aspect_column
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
  
  #----- Input Files
  proteins_tbl_orig <- vroom::vroom(proteins_file, delim="\t")
  
  #----- Compile positive proteins list
  if("UNIPROT_GENENAME" %in% colnames(proteins_tbl_orig)) { 
    positive_proteins <- proteins_tbl_orig %>%
      dplyr::filter(!!rlang::sym(fdr_column_name) < protein_p_val_thresh & 
                      !!rlang::sym(log_fc_column_name) > 0) %>%
      mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, ~str_split(., ":") %>% 
                                                  map_chr(1))) %>% 
      mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
      mutate(gene_name_first = purrr::map_chr(UNIPROT_GENENAME, ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
      mutate(protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, ~str_split(., ":") %>% 
                                                   map_chr(1))) %>%
      group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
      summarise( max_norm_logFC = max(!!rlang::sym(log_fc_column_name))) %>%
      ungroup() %>%
      arrange(comparison, desc(max_norm_logFC))
  } else {
    positive_proteins <- proteins_tbl_orig %>%
      dplyr::filter(!!rlang::sym(fdr_column_name) < protein_p_val_thresh & 
                      !!rlang::sym(log_fc_column_name) > 0) %>%
      mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, ~str_split(., ":") %>% 
                                                  map_chr(1))) %>%
      mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
      group_by(comparison, uniprot_acc_first) %>%
      summarise(max_norm_logFC = max(!!rlang::sym(log_fc_column_name))) %>%
      ungroup() %>%
      arrange(comparison, desc(max_norm_logFC))
  }
  
  vroom::vroom_write(positive_proteins,
                     file.path(output_dir,
                               "all_proteins_with_positive_logFC.tab"),
                     col_names=FALSE)
  
  list_of_comparisons <- positive_proteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) { 
    positive_proteins %>% 
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "all_proteins_with_positive_logFC.tab" ), 
                         col_names = FALSE)})
  
  #----- Compile negative proteins list
  if("UNIPROT_GENENAME" %in% colnames(proteins_tbl_orig)) {
    negative_proteins <- proteins_tbl_orig %>%
      dplyr::filter(!!rlang::sym(fdr_column_name) < protein_p_val_thresh & 
                      !!rlang::sym(log_fc_column_name) < 0) %>%
      mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, ~str_split(., ":") %>% 
                                                  map_chr(1))) %>%
      mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
      mutate(gene_name_first = purrr::map_chr(UNIPROT_GENENAME, ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
      mutate(protein_name_first = purrr::map_chr(`PROTEIN-NAMES`, ~str_split(., ":") %>% 
                                                   map_chr(1))) %>%
      group_by(comparison, uniprot_acc_first, gene_name_first, protein_name_first) %>%
      summarise(min_norm_logFC = min(!!rlang::sym(log_fc_column_name))) %>%
      ungroup() %>%
      arrange(comparison, min_norm_logFC)
  } else {
    negative_proteins <- proteins_tbl_orig %>%
      dplyr::filter(!!rlang::sym(fdr_column_name) < protein_p_val_thresh & 
                      !!rlang::sym(log_fc_column_name) < 0) %>%
      mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, ~str_split(., ":") %>% 
                                                  map_chr(1))) %>%
      mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
      group_by(comparison, uniprot_acc_first) %>%
      summarise(min_norm_logFC = min(!!rlang::sym(log_fc_column_name))) %>%
      ungroup() %>%
      arrange(comparison, min_norm_logFC)
  }
  
  vroom::vroom_write(negative_proteins,
                     file.path(output_dir,
                               "all_proteins_with_negative_logFC.tab"),
                     col_names = FALSE)
  
  list_of_comparisons <- negative_proteins %>% 
    distinct(comparison) %>% 
    pull(comparison)
  
  purrr::walk(list_of_comparisons, ~createDirectoryIfNotExists(file.path(output_dir,. )))
  
  purrr::walk(list_of_comparisons, function(input_comparison) {
    negative_proteins %>%
      dplyr::filter(comparison == input_comparison) %>%
      dplyr::select(uniprot_acc_first) %>%
      vroom::vroom_write(file.path(output_dir,
                                   input_comparison,
                                   "all_proteins_with_negative_logFC.tab" ),
                         col_names = FALSE)})
  
  #----- Compile background proteins list
  background_proteins <- proteins_tbl_orig %>%
    distinct(uniprot_acc) %>%
    mutate(uniprot_acc_first = purrr::map_chr(uniprot_acc, ~str_split(., ":") %>% 
                                                map_chr(1))) %>%
    mutate(uniprot_acc_first = str_replace_all(uniprot_acc_first, "-\\d+$", "")) %>% 
    dplyr::distinct(uniprot_acc_first)
  
  vroom::vroom_write(background_proteins,
                     file.path(output_dir, "background_proteins.tab"),
                     col_names = FALSE)
  
  #----- Compile annotation ID to annotation term name dictionary
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
  
  dictionary <- vroom::vroom(dictionary_file)
  
  id_to_annotation_dictionary <- 
    buildAnnotationIdToAnnotationNameDictionary(input_table=dictionary, 
                                                annotation_column = !!rlang::sym(annotation_column), 
                                                annotation_id_column = !!rlang::sym(annotation_id))
  
  #----- Preparing the enrichment test
  go_annot <- vroom::vroom(annotation_file)
  background_list <- background_proteins
  
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
  
  list_of_comparisons <- negative_proteins %>% 
    bind_rows(positive_proteins) %>%
    distinct(comparison) %>%
    pull(comparison)
  
  #- Tidy up GO aspect list, marked as null if not using GO terms
  go_aspect_list <- NA
  if(!is.null(aspect_column)) { 
    go_aspect_list <- go_annot %>% 
      dplyr::filter(!is.na(!!rlang::sym(aspect_column))) %>%
      distinct(!!rlang::sym(aspect_column)) %>%
      pull(!!rlang::sym(aspect_column)) 
  } else { 
    go_aspect_list <- NA
  }
  
  list_of_genes_list <- list(negative_list = negative_proteins,
                             positive_list = positive_proteins)
  
  input_params <- expand_grid(names_of_genes_list = names( list_of_genes_list), 
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
  
  runOneGoEnrichmentInOutFunctionPartial <- 
    purrr::partial(runOneGoEnrichmentInOutFunction,
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
    mutate(enrichment_results = purrr::pmap(
      list(names_of_genes_list, 
           input_table, 
           go_aspect, 
           input_comparison, 
           min_size, 
           max_size), 
      \(names_of_genes_list, 
        input_table, 
        go_aspect, 
        input_comparison, 
        min_size, 
        max_size) {
        runOneGoEnrichmentInOutFunctionPartial(names_of_genes_list = names_of_genes_list,
                                               input_table = input_table,
                                               go_aspect = go_aspect,
                                               input_comparison = input_comparison,
                                               min_gene_set_size = min_size,
                                               max_gene_set_size = max_size)})) %>% 
    dplyr::select(enrichment_results) %>% 
    unnest(cols = c("enrichment_results"))
  
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
     nrow(enrichment_result) == 0 ) {
    warnings("No enriched terms were identified.")
  } else {
    enrichment_result_add_gene_symbol <- NA
    #- Convert Uniprot accession to gene names
    if(file.exists(uniprot_to_gene_symbol_file)) {
      uniprot_tab_delimited_tbl <- vroom::vroom(file.path(uniprot_to_gene_symbol_file))
      uniprot_to_gene_symbol_dict <- 
        getUniprotAccToGeneSymbolDictionary(uniprot_tab_delimited_tbl, 
                                            !!rlang::sym(protein_id_lookup_column), 
                                            !!rlang::sym(gene_symbol_column), 
                                            !!rlang::sym(protein_id)) 
      convertProteinAccToGeneSymbolPartial <- purrr::partial(convertProteinAccToGeneSymbol,
                                                             dictionary = uniprot_to_gene_symbol_dict)
      enrichment_result_add_gene_symbol <- enrichment_result %>%
        mutate(gene_id_list = str_split(geneID, "/")) %>%
        mutate(gene_symbol = purrr::map_chr(gene_id_list,
                                            convertProteinAccToGeneSymbolPartial)) %>% 
        dplyr::select(-gene_id_list)
    } else {
      enrichment_result_add_gene_symbol <- enrichment_result
    }
    #- generate the output files
    purrr::walk(list_of_comparisons, function(input_comparison) {
      output_file <- paste0(annotation_type, "_table_", input_comparison, ".tab")
      vroom::vroom_write(enrichment_result_add_gene_symbol %>% 
                           dplyr::filter(comparison == input_comparison),
                         file = file.path(output_dir,
                                          input_comparison,
                                          output_file))})
  }
  
}