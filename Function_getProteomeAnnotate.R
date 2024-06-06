#---------- Proteome Annotate
#----- Packages
# suppressPackageStartupMessages({
#   p_load(UniProt.ws)
#   p_load(GO.db)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir <- "results/proteomics/annot_proteins"
# #- Directory path for temporary files
# tmp_dir <- "results/proteomics/cache"
# #- The NCBI taxonomy ID of the organism being investigated (e.g. M.musculus=10090, H.sapien=9606)
# taxonomy_id <- 9606
# #- Results table with values in wider format
# input_wide_file <- "results/proteomics/de_proteins/de_proteins_wide.tsv" 
# #- Results table with values in longer format
# input_long_file <- "results/proteomics/de_proteins/de_proteins_long.tsv"
# #- File to link the cleaned list of accessions to the original list of protein groups 
# #- from MaxQuant file. Also, contains the MaxQuant output row ID
# ids_file <- "results/proteomics/clean_proteins/cleaned_accession_to_protein_group.tab" 
# #- Input file with the protein abundance data
# raw_counts_file <- "data/rawdata/proteinGroups.txt"
# #- Results table with values in wider format
# output_wide_file <- "de_proteins_wide_annot.tsv"
# #- Results table with values in longer format
# output_long_file <- "de_proteins_long_annot.tsv"
# #- Name of the reactome data (Download and save if it does not exists)
# reactome_file <- "UniProt2Reactome.txt"
# #- Name of the uniprot data (Download and save if it does not exists)
# uniprot_file <- "uniprot_data.RDS"

#----- Function
getProteomeAnnotate <- function(
    tmp_dir,
    output_dir,
    taxonomy_id,
    input_wide_file,
    input_long_file,
    ids_file,
    raw_counts_file,
    output_wide_file,
    output_long_file,
    reactome_file,
    uniprot_file
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
  #- Read file with results table with values in longer format
  de_proteins_longer <- vroom::vroom(input_long_file)
  #- Read file to link the cleaned list of accessions to the original list of protein groups from MaxQuant
  ids_tbl <- vroom::vroom(ids_file)
  #- Read file with the protein abundance data
  dat_tbl <- vroom::vroom(raw_counts_file)
  dat_cln <-  janitor::clean_names(dat_tbl) %>%
    dplyr::rename(gene_names_maxquant = "gene_names")
  colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids" )
  
  #----- reactome_file
  reactome_file <- file.path(tmp_dir, reactome_file)
  #- Download Reactome UniProt to pathways file
  if(!file.exists(reactome_file)){
    download.file(url="https://reactome.org/download/current/UniProt2Reactome.txt", 
                  destfile=reactome_file)
  }
  #- Reading Reactome UniProt to pathways file
  reactome_map <- vroom::vroom(reactome_file ,
                               col_names = c("uniprot_acc", 
                                             "reactome_id", 
                                             "url", 
                                             "reactome_term", 
                                             "evidence", 
                                             "organism"))
  
  #----- Get the best UniProt accession per row
  cleanIsoformNumber <- function(string ) {
    # "Q8K4R4-2"
    str_replace(string, "-\\d+$", "")
  }
  uniprot_acc_tbl <- de_proteins_longer %>%
    mutate(uniprot_acc_copy = uniprot_acc) %>%
    separate_rows(uniprot_acc_copy, sep=":") %>%
    mutate(join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy)) %>%
    dplyr::distinct(uniprot_acc, join_uniprot_acc) %>%
    group_by(uniprot_acc) %>%
    mutate(acc_order_id = row_number()) %>%
    ungroup()
  
  #----- Download information from UniProt
  #- The UniProt.ws::select function limits the number of keys queried to 100. 
  #- This gives a batch number for it to be queried in batches
  batchQueryEvidenceHelper <- function(uniprot_acc_tbl, 
                                       uniprot_acc_column){
    #- 100 is the maximum number of queries at one time
    all_uniprot_acc <- uniprot_acc_tbl |>
      dplyr::select({ { uniprot_acc_column } }) |>
      mutate(Proteins = str_split({ { uniprot_acc_column } }, ";")) |>
      unnest(Proteins) |>
      distinct() |>
      arrange(Proteins) |>
      mutate(Proteins = cleanIsoformNumber(Proteins)) |>
      dplyr::mutate(round = ceiling(row_number() / 100))  
  }
  
  # Filter for a batch and run analysis on that batch of uniprot accession keys only.
  subsetQuery <- function(data, 
                          subset, 
                          accessions_col_name, 
                          uniprot_handle, 
                          uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", 
                                              "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                          uniprot_keytype = "UNIPROTKB") {
    # print(subset)
    my_keys <- data |>
      dplyr::filter(round == subset) |>
      pull({ { accessions_col_name } })
    # print(my_keys)
    # print(uniprot_keytype)
    UniProt.ws::select(up,
                       keys = my_keys,
                       columns = uniprot_columns,
                       keytype = uniprot_keytype)
  }
  
  #- Run evidence collection online, giving a table of keys (uniprot_acc_tbl) and 
  #- the column name (uniprot_acc_column)
  batchQueryEvidence <- function(uniprot_acc_tbl, 
                                 uniprot_acc_column, 
                                 uniprot_handle,
                                 uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", 
                                                     "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                                 uniprot_keytype = "UNIPROTKB"){
    # uniprot_evidence_levels <- c("Evidence at protein level",
    #                              "Evidence at transcript level",
    #                              "Inferred from homology",
    #                              "Predicted",
    #                              "Uncertain",
    #                              NA)
    all_uniprot_acc <- batchQueryEvidenceHelper(uniprot_acc_tbl,
                                                {{uniprot_acc_column}})
    
    partial_subset_query <- partial(subsetQuery,
                                    data = all_uniprot_acc,
                                    accessions_col_name = {{uniprot_acc_column}},
                                    uniprot_handle = uniprot_handle,
                                    uniprot_columns = uniprot_columns,
                                    uniprot_keytype = uniprot_keytype)
    
    rounds_list <- all_uniprot_acc |>
      distinct(round) |>
      arrange(round) |>
      pull(round)
    
    all_uniprot_evidence <- purrr::map(rounds_list, \(x){ partial_subset_query(subset = x) }) |>
      bind_rows()
    
    return(all_uniprot_evidence)
  }
  
  uniprot_file <- file.path(tmp_dir, uniprot_file)
  
  if(!file.exists(uniprot_file)){
    up <- UniProt.ws(taxId = taxonomy_id)
    list_of_sp_columns <- c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", 
                            "PROTEIN-NAMES", "LENGTH", "ENSEMBL", "GO-ID", 
                            "KEYWORDS", "protein_existence", "annotation_score", 
                            "reviewed", "gene_names", "protein_name", "length", 
                            "xref_ensembl", "go_id", "keyword")
    up_cls <- unlist(columns(up))
    list_intersect <- intersect(list_of_sp_columns,up_cls)
    # if(length(setdiff(list_of_sp_columns,list_intersect)) > 0)
    # {
    #   logwarn("UniProt fields not found: %s",paste(setdiff( list_of_sp_columns,list_intersect),sep=", "))
    # }
    my_keytype <- "UniProtKB"
    if("UNIPROTKB" %in% keytypes(up)){
      my_keytype <- "UNIPROTKB"
    }
    
    uniprot_dat <- batchQueryEvidence(uniprot_acc_tbl, 
                                      join_uniprot_acc, 
                                      uniprot_handle = up,
                                      uniprot_columns = list_intersect, 
                                      uniprot_keytype = my_keytype)
    
    if(my_keytype == "UniProtKB"){
      uniprot_dat <- uniprot_dat %>%
        dplyr::select(-From) %>%
        dplyr::rename(UNIPROTKB = "Entry",
                      EXISTENCE = "Protein.existence",
                      SCORE = "Annotation",
                      REVIEWED = "Reviewed",
                      GENENAME = "Gene.Names",
                      `PROTEIN-NAMES` = "Protein.names",
                      LENGTH = "Length",
                      ENSEMBL = "Ensembl",
                      `GO-ID` = "Gene.Ontology.IDs",
                      KEYWORDS = "Keywords")
    }
    saveRDS(uniprot_dat, uniprot_file)
  }
  uniprot_dat <- readRDS(uniprot_file)
  
  #----- Merge with Gene Ontology terms
  goterms <- Term(GOTERM)
  gotypes <- Ontology(GOTERM)
  
  uniprotGoIdToTerm <- function(uniprot_dat, 
                                sep = "; ", 
                                goterms, 
                                gotypes){
    uniprot_acc_to_go_id <- uniprot_dat |>
      dplyr::distinct(UNIPROTKB, `GO-ID`) |>
      separate_rows(`GO-ID`, sep = sep) |>
      dplyr::distinct(UNIPROTKB, `GO-ID`) |>
      dplyr::filter(!is.na(`GO-ID`))
    
    go_term_temp <- uniprot_acc_to_go_id |>
      dplyr::distinct(`GO-ID`) |>
      mutate(go_term = purrr::map_chr(`GO-ID`, function(x) { 
        if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
      mutate(go_type = purrr::map_chr(`GO-ID`, function(x) { 
        if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
      mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                                 go_type == "CC" ~ "go_cellular_compartment",
                                 go_type == "MF" ~ "go_molecular_function"))
    
    uniprot_acc_to_go_term <- uniprot_acc_to_go_id |>
      left_join(go_term_temp, by = c("GO-ID" = "GO-ID")) |>
      dplyr::filter(!is.na(go_term)) |>
      group_by(UNIPROTKB, go_type) |>
      summarise(go_term = paste(go_term, collapse = "; ")) |>
      ungroup() |>
      pivot_wider(id_cols = "UNIPROTKB",
                  names_from = go_type,
                  values_from = go_term)
    
    output_uniprot_dat <- uniprot_dat |>
      left_join(uniprot_acc_to_go_term, by = c("UNIPROTKB" = "UNIPROTKB")) |>
      relocate(KEYWORDS, .before = "GO-ID")
    
    return(output_uniprot_dat)
  }
  
  uniprot_dat_cln <- uniprotGoIdToTerm(uniprot_dat, 
                                       sep="; ", 
                                       goterms, 
                                       gotypes)
  
  uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
    left_join(uniprot_dat_cln, 
              by = c("join_uniprot_acc" = "UNIPROTKB")) %>%
    arrange(uniprot_acc, acc_order_id) %>%
    group_by(uniprot_acc) %>%
    summarise(across(.cols=setdiff(colnames(uniprot_dat_cln), "UNIPROTKB"), ~paste(., collapse=":"))) %>%
    ungroup() %>%
    dplyr::rename(UNIPROT_GENENAME = "GENENAME")
  
  #----- Add reactome pathways annotation
  reactome_term_tbl <- uniprot_acc_tbl %>%
    left_join(reactome_map, 
              by = c("join_uniprot_acc" = "uniprot_acc")) %>%
    dplyr::filter(reactome_term != "NA") %>%
    group_by(uniprot_acc, join_uniprot_acc) %>%
    summarise(reactome_term = paste(reactome_term, collapse="; ")) %>%
    ungroup() %>%
    mutate(reactome_term = str_replace_all(reactome_term , ":", "-")) %>%
    group_by(uniprot_acc) %>%
    summarise(reactome_term = paste(reactome_term, collapse=":")) %>%
    ungroup()
  
  #----- Output longer format results table with protein annotation
  de_proteins_longer_annot <- de_proteins_longer %>%
    dplyr::mutate(uniprot_acc_first = str_split(uniprot_acc, ":") |> 
                    purrr::map_chr(1)) |> 
    dplyr::relocate(uniprot_acc_first, .before="uniprot_acc") |>
    left_join(ids_tbl, 
              by = c("uniprot_acc" = "uniprot_acc")) %>%
    left_join(uniprot_dat_multiple_acc, 
              by = c("uniprot_acc" = "uniprot_acc")) %>%
    left_join(reactome_term_tbl, 
              by = c("uniprot_acc" = "uniprot_acc")) %>%
    mutate(maxquant_row_id = as.character(maxquant_row_id)) %>%
    left_join(dat_cln %>% 
                mutate(id = as.character(id)), 
              by=c("maxquant_row_id" = "id", "protein_ids" = "protein_ids")) %>%
    arrange(comparison, q.mod, log2FC) %>%
    distinct() %>%
    relocate(UNIPROT_GENENAME, `PROTEIN-NAMES`, .before="left_group" )
  
  vroom::vroom_write(de_proteins_longer_annot, 
                     file.path(output_dir, output_long_file))
  
  list_of_long_columns <- intersect(colnames(de_proteins_longer_annot), 
                                    c("uniprot_acc", "protein_ids", "protein_names",
                                      "ENSEMBL", "PROTEIN-NAMES", "KEYWORDS",
                                      "GO-ID", "go_biological_process", 
                                      "go_cellular_compartment", "go_molecular_function", 
                                      "reactome_term", "majority_protein_ids", 
                                      "fasta_headers", "evidence_ids", "ms_ms_ids"))
  
  writexl::write_xlsx(de_proteins_longer_annot %>% 
                        mutate_at(list_of_long_columns, ~substr(., 1, 32760)),
                      file.path(output_dir, str_replace(output_long_file, "\\..*", ".xlsx")))
}