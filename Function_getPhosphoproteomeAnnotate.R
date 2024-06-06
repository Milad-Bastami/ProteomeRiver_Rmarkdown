#---------- Phosphoproteome Annotate
#----- Packages
# suppressPackageStartupMessages({
#   p_load(UniProt.ws)
#   p_load(GO.db)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/annot_phos"
# #- Directory path for temporary files
# tmp_dir <- "results/phosphoproteomics/cache"
# #- The NCBI taxonomy ID of the organism being investigated (e.g. M.musculus=10090, H.sapien=9606)
# taxonomy_id <- 9606
# #- Results table with values in wider format
# input_wide_file <- "results/phosphoproteomics/de_phos/de_phos_wide.tsv" 
# #- Results table with values in longer format
# input_long_file <- "results/phosphoproteomics/de_phos/de_phos_long.tsv"
# #- Input file with the protein abundance data
# raw_counts_file <- "results/phosphoproteomics/clean_phos/sum_phoshpsites.tsv"
# #- Results table with values in wider format
# output_wide_file <- "de_phos_wide_annot.tsv"
# #- Results table with values in longer format
# output_long_file <- "de_phos_long_annot.tsv"
# #- File path for location of PhosphoSitePlus DB data files (https://www.phosphosite.org/staticDownloads)
# phosphosite_db_dir <- "data/HomoSapiens/phosphosites"
# #- Find other PTM near within +/- this many number of residues from the phosphosites
# near_ptm_num_residues <- 5
# #- Name of the reactome data (Download and save if it does not exists)
# reactome_file <- "reactome_data.txt"
# #- Name of the uniprot data (Download and save if it does not exists)
# uniprot_file <- "uniprot_data_phos.RDS"

#----- Function
getPhosphoproteomeAnnotate <- function(
    tmp_dir,
    output_dir,
    taxonomy_id,
    input_wide_file,
    input_long_file,
    phosphosite_db_dir,
    raw_counts_file,
    output_wide_file,
    output_long_file,
    near_ptm_num_residues,
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
  ks_file <- file.path(phosphosite_db_dir, 
                       "Kinase_Substrate_Dataset")
  reg_sites_file <- file.path(phosphosite_db_dir, 
                              "Regulatory_sites" )
  disease_file <- file.path(phosphosite_db_dir, 
                            "Disease-associated_sites")
  
  list_of_ptm_file_names <- c("Acetylation_site_dataset",
                              "Methylation_site_dataset",
                              "O-GalNAc_site_dataset",
                              "O-GlcNAc_site_dataset",
                              "Phosphorylation_site_dataset",
                              "Sumoylation_site_dataset",
                              "Ubiquitination_site_dataset")
  
  list_of_ptm_files <- purrr::map_chr(list_of_ptm_file_names, 
                                      ~file.path(phosphosite_db_dir, .))
  
  #----- Read abundance table
  abundance_tbl <- vroom::vroom(raw_counts_file)
  
  #----- Read differentially abundant phosphopeptides table in long format
  de_phos_long <- vroom::vroom(input_long_file) |>
    mutate(sites_id_copy = sites_id, .after = "sites_id") |>
    separate(sites_id_copy, sep="!", into = c("uniprot_acc", "gene_name", "position", "sequence")) |>
    mutate(residue= purrr::map_chr(sequence, 
                                   ~{ str_replace_all( ., "[A-Z_]{7}(.)[A-Z_]{7}([\\:;\\|]*)", "\\1\\2") })) |>
    dplyr::relocate(residue, .before="position")
  
  #----- Read PhosphoSitePlus (PSP) kinase-substrate table
  ks_tbl <- vroom::vroom(ks_file, skip=3) |>
    mutate(SUB_MOD_RSD_CLN = str_replace_all(SUB_MOD_RSD, "([A-Z])(\\d+)", "\\1 \\2")) |>
    separate(SUB_MOD_RSD_CLN, into = c("residue", "position"))
  
  #----- Read PSP regulatory sites table
  reg_sites_tbl <- vroom::vroom(reg_sites_file, skip=3) |>
    mutate(MOD_RSD_CLN = str_replace_all( MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3")) |>
    separate(MOD_RSD_CLN, sep=" ", into=c("residue", "position", "ptm_type")) |>
    relocate(ACC_ID, residue, position, ptm_type, .before="GENE") |>
    dplyr::select(-`...21`) |>
    dplyr::rename(REG_SITES_PMIDs = "PMIDs") |>
    dplyr::rename(REG_SITES_NOTES = "NOTES")
  
  #----- Read PSP disease table
  disease_tbl <- vroom::vroom(disease_file, skip=3) |>
    dplyr::select(ACC_ID, MOD_RSD, DISEASE, ALTERATION, NOTES) |>
    mutate(MOD_RSD_CLN = str_replace_all( MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3")) |>
    separate(MOD_RSD_CLN, sep=" ", into = c("residue", "position", "ptm_type")) |>
    dplyr::rename( DISEASE_NOTES = "NOTES") |>
    dplyr::select(-MOD_RSD, -ptm_type)
  
  #----- Read PSP post-translational modification (PTM) tables
  ptm_tbl <- vroom::vroom(list_of_ptm_files, skip=3, id="ptm") |>
    dplyr::select(ACC_ID, MOD_RSD, ptm) |>
    mutate(MOD_RSD_CLN = str_replace_all(MOD_RSD, "([A-Z])(\\d+)-(.*)", "\\1 \\2 \\3")) |>
    separate(MOD_RSD_CLN, sep=" ", into = c("residue", "position", "ptm_type")) |>
    mutate(uniprot_acc = ACC_ID) |>
    dplyr::mutate(position = as.integer(position)) |>
    mutate(ptm = purrr::map_chr(ptm, ~{ temp_vec <- str_split(., "/")[[1]]
    temp_vec[length(temp_vec)] })) |>
    dplyr::mutate( ptm = str_replace_all(ptm, "_site_dataset", ""))
  
  #----- Find other PTM nearby +/- %s amino acid wihtin the phosphorylation sites
  get_nearby_ptm <- de_phos_long |>
    dplyr::distinct(sites_id, uniprot_acc,  position) |>
    separate_rows(uniprot_acc, position, sep="\\:") |>
    separate_rows(uniprot_acc, position, sep="\\|") |>
    separate_rows(uniprot_acc, position, sep=";") |>
    mutate(position = str_replace_all(position, "\\(|\\)", "") |> 
             purrr::map_int(as.integer)) |>
    mutate(residue_window = purrr::map(position, 
                                       ~seq(from = as.integer( .) - near_ptm_num_residues, 
                                            to = as.integer( .) + near_ptm_num_residues, 
                                            by = 1))) |>
    unnest(residue_window) |>
    left_join(ptm_tbl |> 
                mutate(ptm_count = 1) |>
                dplyr::select(-residue), 
              by=c( "uniprot_acc" = "uniprot_acc", "residue_window" = "position")) |>
    distinct() |>
    dplyr::filter(!(residue_window == position & ptm_type == "p") ) |> 
    ## Avoid counting the same phosphorylation site as the query itself
    dplyr::select(-MOD_RSD, -uniprot_acc, -position, -ptm_type)
  
  nearby_ptm_count <- get_nearby_ptm |>
    dplyr::select(-ACC_ID) |>
    dplyr::filter(!is.na(ptm_count)) |>
    group_by(across(.cols = setdiff(colnames(get_nearby_ptm), 
                                    c("ptm_count", "residue_window", "ACC_ID")))) |>
    distinct() |>
    summarise(ptm_count = sum(ptm_count)) |>
    ungroup() |>
    distinct() |>
    pivot_wider(id_cols = sites_id,
                values_from = ptm_count,
                names_from = "ptm",
                names_prefix = paste0("nearby_+/-", near_ptm_num_residues, "_"))
  
  #----- Number of Phos Sites
  num_phos_sites <- de_phos_long |>
    dplyr::distinct(sites_id, position) |>
    dplyr::mutate(num_sites = str_split(position, ":") |>
                    purrr::map_chr(1) |>
                    str_split( "\\|") |> 
                    purrr::map_chr(1) |>
                    str_split(";") |>
                    purrr::map_int(length)) |>
    dplyr::select(-position)
  
  #----- Phosphosite Table
  phosphosite_plus_tbl <- de_phos_long |>
    dplyr::distinct(sites_id, uniprot_acc, residue, position, sequence) |>
    separate_rows(uniprot_acc, position, sequence, residue, sep="\\:") |>
    separate_rows(uniprot_acc, position, sequence, residue, sep="\\|") |>
    separate_rows(uniprot_acc, position, sequence, residue, sep=";") |>
    mutate(position = str_replace_all(position, "\\(|\\)", "") |> 
             purrr::map_int(as.integer)) |>
    left_join(reg_sites_tbl |>
                dplyr::filter(ptm_type == "p") |>
                dplyr::mutate( position = purrr::map_int(position, as.integer)),
              by=c("uniprot_acc" = "ACC_ID",
                   "residue" = "residue",
                   "position" = "position")) |>
    left_join(ks_tbl |>
                dplyr::rename(KINASE_GENE = "GENE") |>
                dplyr::select(-DOMAIN, - `SITE_+/-7_AA`) |>
                dplyr::mutate(position = purrr::map_int(position, as.integer)),
              by=c("uniprot_acc" = "SUB_ACC_ID",
                   "residue"="residue",
                   "position" = "position",
                   "SITE_GRP_ID" = "SITE_GRP_ID")) |>
    left_join(disease_tbl |>
                dplyr::mutate( position = purrr::map_int(position, as.integer)),
              by=c("uniprot_acc" = "ACC_ID",
                   "residue" = "residue",
                   "position" = "position")) |>
    group_by(sites_id, uniprot_acc) |>
    summarise(across(.cols = everything(), ~paste(unique(.), collapse="//"))) |>
    ungroup() |>
    group_by(sites_id) |>
    summarise(across(.cols = everything(), ~paste(unique(.), collapse=":"))) |>
    ungroup()
  
  #----- Reactome
  reactome_file <- file.path(tmp_dir, reactome_file)
  
  if(!file.exists(reactome_file)) {
    status <- download.file(url = "https://reactome.org/download/current/UniProt2Reactome.txt", 
                            destfile = reactome_file)
  }
  
  reactome_map <- vroom::vroom(reactome_file ,
                               col_names = c("uniprot_acc", "reactome_id", "url", 
                                             "reactome_term", "evidence", "organism"))
  
  #----- Get the best UniProt accession per row
  cleanIsoformNumber <- function(string ) {
    # "Q8K4R4-2"
    str_replace(string, "-\\d+$", "")
  }
  uniprot_acc_tbl <- de_phos_long |>
    mutate(uniprot_acc_copy = uniprot_acc) |>
    separate_rows(uniprot_acc_copy, sep=":") |>
    mutate(join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy)) |>
    dplyr::distinct(uniprot_acc, join_uniprot_acc) |>
    group_by(uniprot_acc) |>
    mutate(acc_order_id = row_number()) |>
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
  
  #- Download information from UniProt
  uniprot_file <- file.path(tmp_dir, uniprot_file)
  if(!file.exists(uniprot_file)) {
    up <- UniProt.ws(taxId = taxonomy_id )
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
    if("UNIPROTKB" %in% keytypes(up)) {
      my_keytype <- "UNIPROTKB"
    }
    uniprot_dat <- batchQueryEvidence(uniprot_acc_tbl, 
                                      join_uniprot_acc, 
                                      uniprot_handle = up,
                                      uniprot_columns = list_intersect, 
                                      uniprot_keytype = my_keytype)
    if(my_keytype == "UniProtKB") {
      uniprot_dat <- uniprot_dat |>
        dplyr::select(-From) |>
        dplyr::rename(UNIPROTKB = "Entry",
                      EXISTENCE = "Protein.existence",
                      SCORE = "Annotation",
                      REVIEWED = "Reviewed",
                      GENENAME = "Gene.Names",
                      `PROTEIN-NAMES` = "Protein.names",
                      LENGTH = "Length",
                      ENSEMBL = "Ensembl",
                      `GO-ID` = "Gene.Ontology.IDs",
                      KEYWORDS   = "Keywords")
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
  de_phos_long_annot <- de_phos_long |>
    dplyr::mutate(uniprot_acc_first = str_split(uniprot_acc, ":") |>
                    purrr::map_chr(1)) |>
    dplyr::relocate(uniprot_acc_first, .before="uniprot_acc") |>
    left_join(num_phos_sites, 
              by = c("sites_id" = "sites_id")) |>
    left_join(uniprot_dat_multiple_acc, 
              by = c("uniprot_acc" = "uniprot_acc")) |>
    left_join(reactome_term_tbl, 
              by = c("uniprot_acc" = "uniprot_acc")) |>
    left_join(phosphosite_plus_tbl |>
                dplyr::select(-uniprot_acc, -position, -residue, -sequence),
              by = c("sites_id" = "sites_id")) |>
    left_join(abundance_tbl |>
                dplyr::select( sites_id, maxquant_row_ids ),
              by = c("sites_id" = "sites_id")) |>
    relocate(maxquant_row_ids, .after="sites_id") |>
    left_join(nearby_ptm_count,
              by = c("sites_id" = "sites_id")) |>
    arrange(comparison, q.mod, log2FC) |>
    distinct()
  
  vroom::vroom_write(de_phos_long_annot, 
                     file.path(output_dir, output_long_file))
  
  list_of_long_columns <- intersect(colnames(de_phos_long_annot), 
                                    c("protein_names", "ENSEMBL", "PROTEIN-NAMES", 
                                      "KEYWORDS", "GO-ID", "go_biological_process",
                                      "go_cellular_compartment", "go_molecular_function", 
                                      "reactome_term", "majority_protein_ids"))
  
  writexl::write_xlsx(de_phos_long_annot |>
                        mutate_at(list_of_long_columns, ~substr(., 1, 32760)),
                      file.path(output_dir,  
                                str_replace(output_long_file, "\\..*", ".xlsx")))
}