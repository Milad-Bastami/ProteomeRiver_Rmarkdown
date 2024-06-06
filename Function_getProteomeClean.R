#----- Function
getProteomeClean <- function(
    tmp_dir,
    output_dir,
    output_counts_file,
    accession_record_file,
    fasta_meta_file,
    fasta_file,
    raw_counts_file,
    column_pattern,
    group_pattern,
    pattern_suffix,
    extract_patt_suffix,
    razor_unique_peptides_group_thresh,
    unique_peptides_group_thresh,
    remove_more_peptides
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
  
  #----- Load Raw Counts File
  dat_tbl <- vroom::vroom(raw_counts_file)
  
  #----- Clean Raw Counts Header
  dat_cln <- janitor::clean_names(dat_tbl)
  colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids")
  
  #----- Prepare regular expressions to select counts data columns
  #- pattern_suffix <- "_\\d+"
  #- extract_patt_suffix <- "_(\\d+)"
  if (group_pattern != "") {
    pattern_suffix <- paste(pattern_suffix, tolower(group_pattern), sep = "_")
    extract_patt_suffix <- paste0(extract_patt_suffix, "_(", tolower(group_pattern), ")")
  }
  
  #- column_pattern <- "reporter_intensity_corrected_\\d+"
  column_pattern <- tolower(paste0(make_clean_names(column_pattern), pattern_suffix))
  
  #- extract_replicate_group <- "reporter_intensity_corrected_d_(\\d+)"
  extract_replicate_group <- tolower(paste0(make_clean_names(column_pattern), 
                                            extract_patt_suffix))
  
  #----- Prepare column names associated with peptide counts
  #- razor_unique_peptides_group_col <- "razor_unique_peptides"
  #- unique_peptides_group_col <- "unique_peptides"
  razor_unique_peptides_group_col <- "razor_unique_peptides"
  unique_peptides_group_col <- "unique_peptides"
  if (group_pattern != "") {
    razor_unique_peptides_group_col <- paste0("razor_unique_peptides_", tolower(group_pattern))
    unique_peptides_group_col <- paste0("unique_peptides_", tolower(group_pattern))
  }
  
  #----- Reading the FASTA file and saving the meta-data file
  #- Extract Fasta Feilds
  getFastaFields <- function(string, pattern) {
    field_found <- str_detect({{string}}, paste0(pattern, "="))
    extract_data <- NA_character_
    if( field_found ) {
      extract_data <- str_replace_all({{string}}, 
                                      paste0("(.*)", pattern, "=(.*?)(\\s..=.*|$)"), "\\2")
    }
    case_when(field_found ~ extract_data,
              TRUE ~ NA_character_  )
  }
  
  #- Clean Isoform Number
  #- clean_isoform_number("Q8K4R4-2")
  cleanIsoformNumber <- function(string) {
    str_replace(string, "-\\d+$", "")
  }
  
  #- Parse FASTA object from seqinr
  parseFastaObject <- function(aa_seq) {
    #- Convert fasta list to data.frame
    #- db, uniprot_acc, uniprot_id, species, tax_id, gene_name, protein_evidence, sequence_version
    accession_tab <-  data.frame(header = names(aa_seq)) %>%
      separate(header, into=c("db", "uniprot_acc", "description"), sep="\\|") %>%
      mutate(uniprot_id = str_replace( description, "(.*?)\\s(.*)", "\\1" ) ) %>%
      mutate(OS = purrr::map_chr(description, ~getFastaFields(., "OS")))  %>%
      mutate(OX = purrr::map_int(description, ~as.integer(getFastaFields(., "OX")))) %>%
      mutate(GN = purrr::map_chr(description, ~getFastaFields(., "GN"))) %>%
      mutate(GN = ifelse( is.na(GN), "", GN)) %>%
      mutate(PE = purrr::map_int(description, ~as.integer(getFastaFields(., "PE")))) %>%
      mutate(SV = purrr::map_int(description, ~as.integer(getFastaFields(., "SV")))) %>%
      dplyr::select(-description) %>%
      dplyr::rename(species = "OS",
                    tax_id = "OX",
                    gene_name = "GN",
                    protein_evidence = "PE",
                    sequence_version = "SV")
    #- Add columns, is_isoform, isoform_num, cleaned_acc, status  
    acc_detail_tab <- accession_tab %>%
      mutate(is_isoform = case_when(str_detect(uniprot_acc, "-\\d+") ~ "Isoform",
                                    TRUE ~ "Canonical")) %>%
      mutate(isoform_num = case_when(is_isoform == "Isoform" ~ str_replace_all(uniprot_acc, 
                                                                               "(.*)(-)(\\d{1,})",
                                                                               "\\3") %>% 
                                       as.numeric, 
                                     is_isoform == "Canonical" ~ 0,
                                     TRUE ~ NA_real_)) %>%
      mutate(cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
      mutate(protein_evidence  = factor(protein_evidence, levels =1:5 )) %>%
      mutate(status = factor( db, levels =c( "sp", "tr"), labels=c("reviewed", "unreviewed"))) %>%
      mutate(is_isoform = factor(is_isoform, levels =c("Canonical", "Isoform")))
    #- Return data.frame
    return(acc_detail_tab)
  }
  
  #- Parse the headers of a Uniprot FASTA file and extract the headers and sequences into a data frame
  parseFastaFile <- function(fasta_file) {
    #- List of all proteins from fasta file
    aa_seqinr <- read.fasta(file = fasta_file,
                            seqtype = "AA",
                            whole.header = TRUE,
                            as.string = TRUE)
    #- Convert fasta list to data.frame
    acc_detail_tab <- parseFastaObject(aa_seq = aa_seqinr)
    #- Rename to Uniprot ID
    names(aa_seqinr) <- str_match(names(aa_seqinr), "(sp|tr)\\|(.+?)\\|(.*)\\s+")[,3]
    #- Add seq and seq_length to data.frame 
    aa_seq_tbl <- acc_detail_tab %>%
      mutate(seq = map_chr(aa_seqinr, 1)) %>%
      mutate(seq_length = purrr::map_int(seq, str_length))
    #- Return data.frame
    return(aa_seq_tbl)
  }
  
  #- fasta_meta_file <- "results/proteomics/cache/aa_seq_tbl.RDS"
  fasta_meta_file <- file.path(tmp_dir, fasta_meta_file)
  if (file.exists(fasta_meta_file)) {
    aa_seq_tbl <- readRDS(fasta_meta_file)
  } else {
    aa_seq_tbl <- parseFastaFile(fasta_file)
    saveRDS(aa_seq_tbl, fasta_meta_file)
  }
  
  #----- Filtering counts table 
  #- Add the row id column and create a column containing the cleaned  peptide
  evidence_tbl <- dat_cln %>%
    mutate(maxquant_row_id = id)
  
  #- Filter count by number of peptides 
  #- Remove decoy proteins and protein contaminants
  select_columns <- evidence_tbl %>%
    dplyr::select(maxquant_row_id,
                  protein_ids,
                  !!rlang::sym(razor_unique_peptides_group_col),
                  !!rlang::sym(unique_peptides_group_col),
                  reverse,
                  potential_contaminant,
                  matches(column_pattern))
  
  #- Remove reverse decoy peptides and contaminant peptides
  if (remove_more_peptides == TRUE) {
    remove_reverse_and_contaminant <- select_columns %>% 
      dplyr::filter(is.na(reverse) &
                      is.na(potential_contaminant)) %>%
      dplyr::filter(!str_detect(protein_ids, "^CON__") &
                      !str_detect(protein_ids, "^REV__")) %>%
      dplyr::filter(!str_detect(protein_ids, "CON__") &
                      !str_detect(protein_ids, "REV__"))
  } else {
    remove_reverse_and_contaminant <- select_columns %>% 
      dplyr::filter(is.na(reverse) &
                      is.na(potential_contaminant)) %>%
      dplyr::filter(!str_detect(protein_ids, "^CON__") &
                      !str_detect(protein_ids, "^REV__"))
  }
  
  #- Melt protein_ids column
  helper_unnest_unique_and_razor_peptides <- remove_reverse_and_contaminant %>%
    dplyr::mutate(protein_ids = str_split(protein_ids, ";")) %>%
    dplyr::mutate(!!rlang::sym(razor_unique_peptides_group_col) := 
                    str_split(!!rlang::sym(razor_unique_peptides_group_col), ";")) %>%
    dplyr::mutate(!!rlang::sym(unique_peptides_group_col) := 
                    str_split(!!rlang::sym(unique_peptides_group_col), ";")) %>%
    unnest(cols = c(protein_ids,
                    !!rlang::sym(razor_unique_peptides_group_col),
                    !!rlang::sym(unique_peptides_group_col)))
  
  #- Filter proteins based on thresholds
  evidence_tbl_cleaned <- helper_unnest_unique_and_razor_peptides %>%
    dplyr::filter(!!rlang::sym(razor_unique_peptides_group_col) >= 
                    razor_unique_peptides_group_thresh & 
                    !!rlang::sym(unique_peptides_group_col) >= 
                    unique_peptides_group_thresh)
  
  #----- Identify best UniProt accession per entry
  #- Identify best UniProt accession per entry, 
  #- extract sample number and simplify column header
  chooseBestProteinAccession <- function(input_tbl, 
                                         acc_detail_tab, 
                                         accessions_column, 
                                         row_id_column = uniprot_acc, 
                                         group_id){
    #- 
    join_condition <- rlang::set_names(c(as_name(enquo(row_id_column)), "cleaned_acc"),
                                       c(as_name(enquo(row_id_column)), "cleaned_acc"))
    
    resolve_acc_helper <- input_tbl |>
      dplyr::select( { { group_id } }, { { accessions_column } }) |>
      mutate( { { row_id_column } } := str_split({ { accessions_column } }, ";")) |>
      unnest( { { row_id_column } }) |>
      mutate( cleaned_acc = cleanIsoformNumber({ { row_id_column } })) |>
      left_join( acc_detail_tab,
                 by = join_condition) |>
      dplyr::select( { { group_id } }, one_of(c(as_name(enquo(row_id_column)), "gene_name", "cleaned_acc",
                                                "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"))) |>
      distinct() |>
      arrange( { { group_id } }, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)
    
    score_isoforms <- resolve_acc_helper |>
      mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
      group_by({ { group_id } }, gene_name) |>
      arrange( { { group_id } }, protein_evidence,
               status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) |>
      mutate( ranking = row_number()) |>
      ungroup()
    
    ## For each gene name find the uniprot_acc with the lowest ranking
    my_group_id <- enquo(group_id)
    
    join_names <- rlang::set_names(c(as_name(my_group_id), "ranking", "gene_name"),
                                   c(as_name(my_group_id), "ranking", "gene_name"))
    
    group_gene_names_and_uniprot_accs <- score_isoforms |>
      distinct( { { group_id } }, gene_name, ranking) |>
      dplyr::filter(ranking == 1) |>
      left_join(score_isoforms |>
                  dplyr::select({ { group_id } }, ranking, gene_name, uniprot_acc),
                by = join_names) |>
      dplyr::select(-ranking) |>
      group_by({ { group_id } }) |>
      summarise(num_gene_names = n(),
                gene_names = paste(gene_name, collapse = ":"),
                uniprot_acc = paste(uniprot_acc, collapse = ":")) |>
      ungroup() |>
      mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                   TRUE ~ "Multimapped"))
    return(group_gene_names_and_uniprot_accs)
  }
  #- UniProt Accession with Gene List
  accession_gene_name_tbl <- chooseBestProteinAccession(input_tbl = evidence_tbl_cleaned,
                                                        acc_detail_tab = aa_seq_tbl,
                                                        accessions_column = protein_ids,
                                                        row_id_column = uniprot_acc,
                                                        group_id = maxquant_row_id)
  #- UniProt Accession with Gene List and Protein IDs
  accession_gene_name_tbl_record <- accession_gene_name_tbl %>%
    left_join(evidence_tbl %>% 
                dplyr::select(maxquant_row_id, protein_ids), by = c("maxquant_row_id"))
  
  #- UniProt Accession with sample expressions
  evidence_tbl_filt <- evidence_tbl_cleaned |>
    inner_join(accession_gene_name_tbl |>
                 dplyr::select(maxquant_row_id, uniprot_acc), by = "maxquant_row_id") |>
    dplyr::select(uniprot_acc, matches(column_pattern), -contains(c("razor", "unique"))) |>
    distinct()
  
  #- Extraction Pattern
  if (group_pattern != "") {
    extraction_pattern <- "\\1_\\2"
  } else {
    extraction_pattern <- "\\1"
  }
  
  colnames(evidence_tbl_filt) <- 
    str_replace_all(colnames(evidence_tbl_filt), 
                    tolower(extract_replicate_group), extraction_pattern) %>%
    toupper( ) %>%
    str_replace_all("UNIPROT_ACC", "uniprot_acc")
  
  #- Fix Column names
  colnames(evidence_tbl_filt) <- gsub(pattern = "REPORTER_INTENSITY_CORRECTED_", 
                                      replacement = "", 
                                      x = colnames(evidence_tbl_filt))
  
  #----- Save Data
  vroom::vroom_write(evidence_tbl_filt, 
                     file.path(output_dir, output_counts_file))
  
  vroom::vroom_write(accession_gene_name_tbl_record, 
                     file.path(output_dir, accession_record_file))
  
  vroom::vroom_write(data.frame(sample_names = colnames(evidence_tbl_filt)[-1]), 
                     file.path(output_dir, "sample_names.tab"))
  
  # #- Array to store number of proteins in each step
  # num_proteins_remaining <- rep(NA, 3)
  # names(num_proteins_remaining) <- 
  #   c("Number of proteins in raw unfiltered file", 
  #     "Number of proteins after removing reverse decoy and contaminant proteins", 
  #     paste0("Number of proteins after removing proteins with no. of razor + unique peptides < ",  
  #            razor_unique_peptides_group_thresh, 
  #            " and ",  "no. of unique peptides < ", 
  #            unique_peptides_group_thresh))
  # # Record the number of proteins in raw unfiltered file
  # num_proteins_remaining[1] <- nrow(select_columns)
  # # Record the number of proteins after removing reverse decoy and contaminant proteins
  # num_proteins_remaining[2] <- nrow(remove_reverse_and_contaminant)
  # # Record the number of proteins after filtering razor + unique peptides 
  # num_proteins_remaining[3] <- nrow(evidence_tbl_filt)
  # # Record the number of proteins remaining after each filtering step into the file 'number_of_proteins_remaining_after_each_filtering_step.tab'
  # num_proteins_remaining_tbl <- data.frame(step=names( num_proteins_remaining), 
  #                                          num_proteins_remaining=num_proteins_remaining)
  # vroom::vroom_write(num_proteins_remaining_tbl, 
  #                    file.path(output_dir, 
  #                              "number_of_proteins_remaining_after_each_filtering_step.tab"))
}