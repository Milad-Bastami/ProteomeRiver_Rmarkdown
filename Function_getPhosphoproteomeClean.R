#---------- Phosphoproteome Clean
#----- Packages
# suppressPackageStartupMessages({
#   p_load(janitor)
#   p_load(tidyverse)
#   p_load(seqinr)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/clean_phos"
# #- Directory path for temporary files
# tmp_dir <- "results/phosphoproteomics/cache"
# #- Input protein sequence FASTA file with UniProt FASTA header format
# fasta_file <- "data/HomoSapiens/fasta/HomoSapiensCanIso_2023_12_06.fasta"
# #- R object storage file that records all the sequence and header information in the FASTA file
# fasta_meta_file <- "aa_seq_tbl.RDS"
# #- Input file with the protein abundance data
# raw_counts_file <- "data/rawdata/evidence.txt"
# #- Numeric value representing the probability threshold for accepting a primary modification site
# site_prob_threshold <- 0.75
# #- Numeric value representing the probability threshold for recovering the abundance data for a 
# #- lower quality modification site, if the same site has already been found as high probability site
# recover_site_prob_thresh <- 0.5
# #- String pattern that matches the abundance data columns. It also facilitate the ID of the sample to be extracted
# col_pattern_string <- "Reporter intensity corrected"
# #- Regular expression pattern to identify columns with abundance values belonging to the experiment
# pattern_suffix <- "_\\d+"
# #- Regular expression pattern to identify columns with abundance values belonging to the experiment
# extract_patt_suffix <- "_(\\d+)"
# #- A string listing all the additional columns to be included (e.g. experiment group column). 
# #- Each column is separated by a comma
# add_cols_string <- ""

#----- Function
getPhosphoproteomeClean <- function(
    tmp_dir,
    output_dir,
    fasta_meta_file,
    fasta_file,
    raw_counts_file,
    col_pattern_string,
    add_cols_string,
    pattern_suffix,
    extract_patt_suffix,
    site_prob_threshold,
    recover_site_prob_thresh
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
  evidence_tbl <- vroom::vroom(raw_counts_file)
  
  #----- Clean Raw Counts Header
  evidence_janitor <- janitor::clean_names(evidence_tbl)
  colnames(evidence_janitor) <- str_replace(colnames(evidence_janitor), "_i_ds", "_ids" )
  
  #----- Column Names
  additional_cols <- str_split(add_cols_string, ",")[[1]]
  col_pattern <- janitor::make_clean_names(col_pattern_string)
  
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
  
  #- fasta_meta_file <- "results/phosphoproteomics/cache/aa_seq_tbl.RDS"
  fasta_meta_file <- file.path(tmp_dir, fasta_meta_file)
  aa_seq_tbl <- NA
  if(!file.exists(file.path(fasta_meta_file))) {
    aa_seq_tbl <- parseFastaFile(fasta_file)
    saveRDS(aa_seq_tbl, file.path(fasta_meta_file))
  } else {
    aa_seq_tbl <- readRDS(file.path(fasta_meta_file))
  }
  
  #----- Add the row id column and create a column containing the cleaned peptide
  addColumnsToEvidenceTbl <- function(evidence_tbl, 
                                      phospho_site_prob_col = phospho_sty_probabilities) {
    evidence_tbl_cleaned <- evidence_tbl %>%
      mutate(evidence_id = (row_number() - 1)) %>%
      # dplyr:::select(one_of(c("evidence_id", evidence_col_to_use %>% pull(Columns)))) %>%
      mutate(cleaned_peptide = str_replace_all({{phospho_site_prob_col}}, "[\\(\\)0-9\\.]", ""))
    return( evidence_tbl_cleaned)
  }
  
  evidence_tbl_cleaned <- addColumnsToEvidenceTbl(evidence_janitor)
  
  #----- Get best accession per entry, work out peptides mapped to multiple genes
  #- Use decision tree to get best accession per phosphosite evidence entry
  chooseBestPhosphositeAccession <- function(input_tbl, 
                                             acc_detail_tab, 
                                             accessions_column, 
                                             group_id) {
    resolve_acc_helper <- input_tbl %>%
      dplyr::select({{group_id}}, {{accessions_column}}, cleaned_peptide) %>%
      mutate(uniprot_acc = str_split({{accessions_column}}, ";")) %>%
      unnest(uniprot_acc)   %>%
      mutate(cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
      left_join(acc_detail_tab,
                by=c("uniprot_acc" = "uniprot_acc",
                     "cleaned_acc" = "cleaned_acc")) %>%
      ## Just a sanity check that the peptide is actually in the sequence
      dplyr::filter(str_detect(seq, cleaned_peptide)) %>%
      dplyr::select({{group_id}}, one_of(c( "uniprot_acc", "gene_name", "cleaned_acc",
                                            "protein_evidence", "status", "is_isoform", 
                                            "isoform_num", "seq_length"))) %>%
      distinct %>%
      arrange({{group_id}}, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)
    
    print(colnames(head(resolve_acc_helper)))
    
    score_isoforms <- resolve_acc_helper %>%
      mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
      group_by({{group_id}}, gene_name) %>%
      arrange({{group_id}}, protein_evidence, status, is_isoform, desc(seq_length), 
              isoform_num, cleaned_acc) %>%
      mutate(ranking = row_number()) %>%
      ungroup()
    
    print( colnames(head(score_isoforms)) )
    
    ## For each gene name find the uniprot_acc with the lowest ranking
    my_group_id <- enquo(group_id)
    
    join_names <-
      rlang::set_names(c(as_name(my_group_id), "ranking", "gene_name"),
                       c(as_name(my_group_id), "ranking", "gene_name"))
    
    group_gene_names_and_uniprot_accs <- score_isoforms %>%
      distinct({{group_id}}, gene_name, ranking) %>%
      dplyr::filter(ranking == 1) %>%
      left_join(score_isoforms %>%
                  dplyr::select({{group_id}}, ranking, gene_name, uniprot_acc),
                by = join_names) %>%
      dplyr::select(-ranking)
    
    # %>%
    #   group_by({{group_id}}) %>%
    #   summarise( num_gene_names = n(),
    #              gene_names = paste( gene_name, collapse=":"),
    #              uniprot_acc = paste( uniprot_acc, collapse=":")) %>%
    #   ungroup() %>%
    #   mutate( is_unique = case_when( num_gene_names == 1 ~ "Unique",
    #                                  TRUE ~ "Multimapped"))
    return( group_gene_names_and_uniprot_accs )
  }
  
  accession_gene_name_tbl <- chooseBestPhosphositeAccession(evidence_tbl_cleaned,
                                                            aa_seq_tbl,
                                                            leading_proteins,
                                                            evidence_id)
  
  #----- Remove peptides without abundance values at all
  removePeptidesWithoutAbundances <- function(evidence_tbl_cleaned, 
                                              col_pattern) {
    sites_to_accept <- evidence_tbl_cleaned %>%
      mutate(across(matches(col_pattern, perl=TRUE), ~.==0)) %>%
      dplyr::filter(!if_all(matches(col_pattern, perl=TRUE), ~. ==TRUE)) %>%
      dplyr::select(evidence_id)
    ## Removing entries where all the "Reporter intensity corrected" rows are zero
    evidence_tbl_filtered <- evidence_tbl_cleaned %>%
      inner_join(sites_to_accept, by = c("evidence_id" = "evidence_id"))
    return(evidence_tbl_filtered)
  }
  
  evidence_tbl_filt <- removePeptidesWithoutAbundances(evidence_tbl_cleaned, col_pattern)
  
  #----- Filter peptides with no intensity across all samples, extract intensity data, extract sites
  #- For all the multi-phosphosites peptide extract their intensity, filter peptide 
  #- with no intensity across all samples, extract site probabilities
  getMaxProb <- function(phosphopeptide, 
                         num_sites=1) {
    pass_thresh <- str_match_all(phosphopeptide, "\\((\\d+\\.*\\d*)\\)") %>% 
      .[[1]]  %>%
      .[,2] %>%
      as.numeric()
    # Try to preserve the order in which the probability is listed in the peptide,
    # while using sort to find the top 'num_sites'
    if(length(pass_thresh) == 0) {
      return(c())
    }
    # Sort from maximum to minimum and then take the first few numbers according to number of sites required
    # Please ensure decreasing is set to TRUE.
    top_site_index <- sort.int(pass_thresh,
                               index.return=TRUE,
                               decreasing=TRUE)$ix[seq_len(num_sites)]
    pass_thresh[sort(top_site_index)] %>%
      keep(~{ !is.na(.) })
  }
  
  getMaxProbFutureMap <- function(phosphopeptide, 
                                  num_sites=1) {
    furrr::future_map2(phosphopeptide, 
                       num_sites,
                       ~{ getMaxProb(.x, .y) })
  }
  
  getBestPosition <- function(phosphopeptide, 
                              num_sites=1) { 
    if(str_detect(phosphopeptide, "p")) {
      stop("Input phosphopetide string should not have little 'p' as characters.")
    }
    
    pass_thresh <- str_match_all(phosphopeptide,
                                 "\\((\\d+\\.*\\d*)\\)") %>%
      .[[1]]  %>%
      .[,2] %>%
      as.numeric()
    
    prob_list <- getMaxProb(phosphopeptide, 
                            num_sites)
    ## I might need to fix this line as if we have two poisition sharing the same maximum score,
    ## we currently only use the first one as best position
    selected_pos <- which( pass_thresh %in% prob_list)
    
    little_p_position <- str_replace_all(phosphopeptide,
                                         "\\(\\d+\\.*\\d*\\)", "p") %>%
      str_locate_all("p") %>%
      .[[1]] %>%
      .[,1]
    
    to_adj_pos <- seq_along(little_p_position)
    clean_pos <- little_p_position - to_adj_pos
    return( clean_pos[selected_pos] )
  }
  
  getBestPositionFutureMap <- function(phosphopeptide, 
                                       num_sites = 1) {
    furrr::future_map2(phosphopeptide, 
                       num_sites,
                       ~{ getBestPosition(.x, .y) })
  }
  
  filterPeptideAndExtractProbabilities <- function(evidence_tbl_cleaned, 
                                                   accession_gene_name_tbl,
                                                   col_pattern = "corrected",
                                                   accession_col = leading_proteins,
                                                   phospho_site_prob_col = phospho_sty_probabilities,
                                                   num_phospho_site_col = phospho_sty) {
    sites_probability_tbl <- evidence_tbl_cleaned %>%
      ## Must be a phosphopeptide (at least one site)
      dplyr::filter({{num_phospho_site_col}} >=1) %>%
      ## Remove REV_ and CON_
      dplyr::filter(!str_detect({{accession_col}}, "^REV__|^CON__")) %>%
      dplyr::mutate(best_phos_prob = getMaxProbFutureMap({{phospho_site_prob_col}},
                                                         {{num_phospho_site_col}})) %>%
      dplyr::filter(map_lgl(best_phos_prob, ~{ length(.) > 0 })) %>%
      dplyr::mutate(best_phos_pos = getBestPositionFutureMap({{phospho_site_prob_col}},
                                                             {{num_phospho_site_col}})) %>%
      ## Avoid cases where there are multiple positions having the same top scores
      dplyr::filter(map2_lgl(best_phos_prob, best_phos_pos, ~{ length(.x) == length(.y) })) %>%
      left_join(accession_gene_name_tbl, by="evidence_id") %>%
      dplyr::select(one_of(c("best_phos_prob", "best_phos_pos", 
                             colnames(accession_gene_name_tbl),
                             colnames(evidence_tbl_cleaned))), 
                    {{phospho_site_prob_col}},
                    {{num_phospho_site_col}})
    
    number_of_rows_without_uniprot_acc <- sites_probability_tbl %>%
      dplyr::filter(is.na(uniprot_acc)) %>%
      nrow()
    
    if(length(number_of_rows_without_uniprot_acc) > 0) {
      warnings(paste("There are", number_of_rows_without_uniprot_acc, 
                     "proteins do not have sequence information in FASTA file. Removing them from analysis"))
    }
    
    sites_probability_filt <- sites_probability_tbl %>%
      dplyr::filter(!is.na(uniprot_acc))
    
    return(sites_probability_filt )
  }
  
  sites_probability_tbl <- filterPeptideAndExtractProbabilities(evidence_tbl_filt,
                                                                accession_gene_name_tbl,
                                                                col_pattern,
                                                                accession_col = leading_proteins,
                                                                phospho_site_prob_col = phospho_sty_probabilities,
                                                                num_phospho_site_col = phospho_sty)
  
  #----- Get the peptide start and end position for each peptide
  addPeptideStartAndEnd <- function(sites_probability_tbl, 
                                    aa_seq_tbl) {
    peptide_start_and_end <- sites_probability_tbl %>%
      left_join(aa_seq_tbl %>% 
                  dplyr::select(uniprot_acc, seq), 
                by = c("uniprot_acc" = "uniprot_acc")) %>%
      mutate(peptide_location = str_locate_all(seq, cleaned_peptide)) %>%
      mutate(pep_start = map_chr(peptide_location, ~paste(.[,"start"], collapse="|"))) %>%
      mutate(pep_end = map_chr(peptide_location, ~paste(.[,"end"], collapse="|")))
    return( peptide_start_and_end)
  }
  
  peptide_start_and_end <- addPeptideStartAndEnd(sites_probability_tbl, 
                                                 aa_seq_tbl)
  
  #----- Get the phosphosites position string
  getPosString <-  function(peptide_start_position, 
                            site_relative_position) {
    a <- peptide_start_position
    b <- site_relative_position
    pos_group <- cross2(a, b) %>%
      map_dbl( ~{sum(unlist(.))-1} )
    pos_mat <- matrix(pos_group,
                      ncol = length(b),
                      nrow = length(a),
                      byrow = FALSE)
    # print( pos_mat)
    if(length(a) > 1) {
      pos_string <- list()
      for(i in seq_len(nrow(pos_mat))) {
        pos_string <- c(pos_string, paste0( "(", paste( pos_mat[i,], collapse=";"), ")"))
      }
      return(paste(pos_string, collapse="|"))
    } else {
      pos_string <- paste(pos_mat[1,], collapse=";")
      return(pos_string)
    }
  }
  
  addPhosphositesPositionsString <- function(peptide_start_and_end) {
    phosphosite_pos_string_tbl <- peptide_start_and_end %>%
      mutate(best_phos_pos_string = map_chr(best_phos_pos, 
                                            ~paste(., collapse=";"))) %>%
      mutate(temp_check_pos = map2(peptide_location, 
                                   best_phos_pos, 
                                   ~{cross2( .x[,"start"] , .y)} )) %>%
      mutate(check_pos = purrr::map(temp_check_pos, 
                                    ~{ map_dbl(., function(x){sum(unlist(x)) -1})} )) %>%
      mutate(protein_site_positions = map2_chr(peptide_location, 
                                               best_phos_pos, 
                                               ~{getPosString(.x[, "start"] , .y) })) %>%
      mutate( best_phos_prob_string = map_chr(best_phos_prob, 
                                              ~paste(., collapse=";")))
    return(phosphosite_pos_string_tbl)
  }
  
  phosphosite_pos_string_tbl <- addPhosphositesPositionsString(peptide_start_and_end)
  
  #----- Get the string listing all the 15-mer sequences, each sequence has the phosphorylation site in the middle
  getXMerString <- function(seq, 
                            uniprot_acc, 
                            position, 
                            padding_length = 7) {
    start <- position - padding_length
    end <- position + padding_length
    seq_end <- str_length(seq)
    # print( seq_end)
    end_padding <- 0
    start_padding <- 0
    if(seq_end < end) {
      end_padding <- position + padding_length - seq_end
      end <- seq_end
    }
    if(start < 1) {
      start_padding <- abs(start - 1)
      start <- 1
    }
    start_padding_string <- paste0(rep("_", start_padding), collapse="")
    end_padding_string <- paste0(rep("_", end_padding), collapse="" )
    my_X_mer_partial <- str_sub(seq, start, end)[[1]]
    my_X_mer <- paste0(start_padding_string,
                       my_X_mer_partial,
                       end_padding_string)
    return(my_X_mer)
  }
  
  getXMersList <-  function(seq, 
                            uniprot_acc,
                            peptide_start_position, 
                            site_relative_position, 
                            padding_length = 7) {
    a <- peptide_start_position
    b <- site_relative_position
    pos_group <- cross2(a, b) %>%
      map_dbl( ~{sum(unlist(.))-1} )
    pos_mat <- matrix(pos_group,
                      ncol = length(b),
                      nrow = length(a),
                      byrow = FALSE)
    # print( uniprot_acc)
    # print( pos_mat)
    # print(as.vector(pos_mat[1,] ) )
    # print( paste( "peptide_start_position = ", paste( peptide_start_position, collapse=";")))
    # print( paste( "site_relative_position = ", paste( site_relative_position, collapse=";")))
    my_Xmers_list <- purrr::map_chr( as.vector(pos_mat[1,] ), 
                                     ~{getXMerString(seq, uniprot_acc, ., padding_length=padding_length)}) %>%
      paste(collapse=";")
    return(my_Xmers_list)
  }
  
  addXMerStrings <- function(phosphosite_pos_string_tbl, 
                             padding_length=7) {
    my_get_X_mers_list <- function(uniprot_acc, 
                                   peptide_location, 
                                   best_phos_pos, 
                                   seq) { 
      getXMersList(seq, 
                   uniprot_acc, 
                   peptide_location, 
                   best_phos_pos, 
                   padding_length = padding_length)
    }
    
    get_15_mer_tbl <- phosphosite_pos_string_tbl %>%
      mutate(phos_15mer_seq = furrr::future_pmap_chr(list(uniprot_acc = uniprot_acc,
                                                          peptide_location = peptide_location,
                                                          best_phos_pos = best_phos_pos,
                                                          seq = seq), 
                                                     my_get_X_mers_list)) %>% 
      distinct()
    return( get_15_mer_tbl)
  }
  
  get_15_mer_tbl <- addXMerStrings(phosphosite_pos_string_tbl, 
                                   padding_length = 7)
  
  #----- Get peptides with at least one phosphosite over threshold. 
  #- Find all peptides with same sites as another peptides that contained at least one phosphosies >= threshold.
  filterByScoreAndGetSimilarPeptides <- function(get_15_mer_tbl, 
                                                 site_prob_threshold, 
                                                 secondary_site_prob_threshold = 0.5, 
                                                 num_phospho_site_col = phospho_sty) {
    #- Find peptide in which at least one phosphosite has one position >= site probability threshold
    all_peptide_and_sites_pass_filter <- get_15_mer_tbl %>%
      dplyr::filter(map2_lgl(best_phos_prob, 
                             {{num_phospho_site_col}},
                             function(x, y) { 
                               return(length(which( x >= site_prob_threshold)) >=  y) })) %>%
      dplyr::select(uniprot_acc, 
                    protein_site_positions) %>%
      distinct() %>%
      dplyr::mutate(protein_site_positions = str_split(protein_site_positions, "[;\\|]")) %>%
      unnest(protein_site_positions) %>%
      dplyr::mutate(protein_site_positions = as.integer(str_replace_all(protein_site_positions, 
                                                                        "\\(|\\)", ""))) %>%
      distinct() %>%
      arrange(uniprot_acc, protein_site_positions) %>%
      group_by(uniprot_acc) %>%
      nest(high_quality_sites = protein_site_positions) %>%
      ungroup() %>%
      mutate(high_quality_sites = purrr::map(high_quality_sites, 
                                             ~{ .$protein_site_positions }))
    #- Find peptide where all major sites probability is > secondary_site_prob_threshold
    peptide_and_pos_pass_filt <- get_15_mer_tbl %>%
      dplyr::filter(map2_lgl(best_phos_prob, 
                             {{num_phospho_site_col}},
                             function(x, y) { return(length(which(x > secondary_site_prob_threshold)) >= y)})) %>%
      dplyr::select(uniprot_acc, protein_site_positions) %>%
      distinct() %>%
      dplyr::mutate(protein_site_positions = str_split(protein_site_positions, 
                                                       "[;\\|]")) %>%
      dplyr::mutate(protein_site_positions = purrr::map(protein_site_positions, 
                                                        ~{ str_replace_all(., "\\(|\\)", "") %>% 
                                                            purrr::map_int(as.integer) })) %>% 
      distinct()
    
    all_filtered_peptide <- peptide_and_pos_pass_filt %>%
      dplyr::inner_join(all_peptide_and_sites_pass_filter, by = "uniprot_acc") %>%
      dplyr::filter(map2_lgl(protein_site_positions,
                             high_quality_sites,
                             function(to_check, hq_sites) { 
                               length(which(to_check %in% hq_sites)) == length(to_check)})) %>%
      dplyr::mutate(protein_site_positions = map_chr(protein_site_positions,
                                                     ~paste(., collapse=";")))
    #- Check that all sites in peptide has been found at least one in set "all_peptide_and_sites_pass_filter"
    #- Keep peptide if there is another peptide that has >0.75 at all of these sites 
    #- and secondary site with highest probability in the same position.
    get_15_mer_tbl_filt <- get_15_mer_tbl %>%
      dplyr::mutate(temp_protein_site_positions = purrr::map_chr(protein_site_positions, 
                                                                 ~{ str_replace_all(., "\\(|\\)", "") %>%
                                                                     str_replace_all( "\\|", ";" )})) %>%
      dplyr::inner_join(all_filtered_peptide, 
                        by = c("uniprot_acc" = "uniprot_acc",
                               "temp_protein_site_positions" = "protein_site_positions")) %>%
      dplyr::select(-temp_protein_site_positions)
    return(get_15_mer_tbl_filt)
  }
  
  get_15_mer_tbl_filt <- filterByScoreAndGetSimilarPeptides(get_15_mer_tbl,
                                                            site_prob_threshold,
                                                            recover_site_prob_thresh)
  
  #----- Pivot the phosphosites to a longer table
  allPhosphositesPivotLonger <- function(get_15_mer_tbl,
                                         additional_cols = c("experiment"),
                                         col_pattern = "Reporter intensity corrected",
                                         pattern_suffix = "_\\d+",
                                         extract_patt_suffix = "_(\\d+)",
                                         phospho_site_prob_col = phospho_sty_probabilities,
                                         num_phospho_site_col = phospho_sty) {
    usual_columns <- c("evidence_id", "uniprot_acc", "gene_name", "sequence", # "gene_names",
                       "protein_site_positions", "phos_15mer_seq", 
                       as_name(enquo(phospho_site_prob_col)), 
                       as_name(enquo(num_phospho_site_col)))
    cols_to_use <- usual_columns
    if(!is.na(additional_cols) & additional_cols != "") {
      cols_to_use <- c(usual_columns, additional_cols)
    }
    all_sites_long <- get_15_mer_tbl %>%
      dplyr::select({{cols_to_use}},
                    matches(paste0(tolower(col_pattern), pattern_suffix), perl = TRUE)) %>%
      pivot_longer(cols = matches(c(paste0(tolower(col_pattern), pattern_suffix)), perl = TRUE),
                   names_to = "replicate",
                   values_to = "value")
    # print(head(all_sites_long))
    # print( col_pattern)
    if(extract_patt_suffix != "") {
      all_sites_long <- all_sites_long %>%
        dplyr::mutate(replicate = str_replace(replicate,
                                              paste0(tolower(col_pattern), extract_patt_suffix),
                                              "\\1") %>% 
                        map_int(as.integer))
    }
    return(all_sites_long)
  }
  
  all_phos_sites_long_tbl <- allPhosphositesPivotLonger(get_15_mer_tbl_filt,
                                                        additional_cols ,
                                                        col_pattern,
                                                        pattern_suffix = pattern_suffix,
                                                        extract_patt_suffix = extract_patt_suffix)
  
  #----- Group peptides from paralog proteins
  groupParalogPeptides <- function(all_sites_long,
                                   additional_cols = c("experiment"),
                                   phospho_site_prob_col = phospho_sty_probabilities,
                                   num_phospho_site_col = phospho_sty) {
    grouping_variables <- c("evidence_id", "replicate", "value", "sequence",
                            as_name(enquo( phospho_site_prob_col)),
                            as_name(enquo( num_phospho_site_col)))
    if(!is.na(additional_cols) & additional_cols != "") {
      grouping_variables <- c(grouping_variables, additional_cols)
    }
    final_select_var <- c("evidence_id", "uniprot_acc", "gene_names", "protein_site_positions", 
                          "phos_15mer_seq", "replicate", "value",
                          as_name(enquo( phospho_site_prob_col)),
                          as_name(enquo( num_phospho_site_col)) )
    if(!is.na(additional_cols) & additional_cols != "") {
      final_select_var <- c(final_select_var,
                            additional_cols )
    }
    #- Group homolog gene names, uniprot_acc, site posiitons
    paralog_sites_long <- all_sites_long %>%
      group_by(across({{grouping_variables}})) %>%
      summarise(gene_names = paste(gene_name, collapse = ":"),
                uniprot_acc = paste(uniprot_acc, collapse=":"),
                protein_site_positions = paste( protein_site_positions, collapse=":"),
                phos_15mer_seq = paste( phos_15mer_seq, collapse=":")) %>%
      ungroup() %>%
      dplyr::select(one_of(final_select_var))
    return(paralog_sites_long)
  }
  
  paralog_sites_long <- groupParalogPeptides(all_phos_sites_long_tbl, 
                                             additional_cols)
  
  #----- Pivot the phosphosites data to a wide format
  allPhosphositesPivotWider <- function(all_phos_sites_long_tbl,
                                        additional_cols = c("experiment"),
                                        phospho_site_prob_col = phospho_sty_probabilities,
                                        num_phospho_site_col = phospho_sty) {
    cols_to_use <- "replicate"
    temp_tbl <- all_phos_sites_long_tbl
    if(!is.na(additional_cols) & additional_cols != "") {
      cols_to_use <-c("replicate", additional_cols)
      temp_tbl <- all_phos_sites_long_tbl %>%
        mutate_at(additional_cols, toupper)
    }
    all_phos_sites_wide_tbl <- temp_tbl %>%
      pivot_wider(id_cols = c(evidence_id, uniprot_acc, gene_names,
                              protein_site_positions, phos_15mer_seq,
                              as_name(enquo( phospho_site_prob_col)),
                              as_name(enquo( num_phospho_site_col))),
                  names_from = all_of(cols_to_use),
                  values_from = value) %>%
      arrange(uniprot_acc, gene_names, protein_site_positions,
              phos_15mer_seq, evidence_id)
    return(all_phos_sites_wide_tbl)
  }
  
  all_phos_sites_wide_tbl <- allPhosphositesPivotWider(paralog_sites_long,
                                                       additional_cols)
  
  #----- Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in long format
  uniquePhosphositesSummariseLongList <- function(all_phos_sites_long_tbl,
                                                  additional_cols = c("experiment")) {
    #- Summarise the input table with a summarisation function
    group_summary <- function( input_tbl, additional_cols, method=mean ) {
      usual_columns <- c( "uniprot_acc",  "gene_names", "protein_site_positions", 
                          "phos_15mer_seq", "replicate")
      cols_to_use <- usual_columns
      if(!is.na(additional_cols) & additional_cols != "") {
        cols_to_use <- c(usual_columns, additional_cols)
      }
      temp_tbl <- input_tbl %>%
        group_by(across({{ cols_to_use }})) %>%
        # uniprot_acc, gene_names, protein_site_positions, phos_15mer_seq, experiment, replicate
        summarise(value = method(value),
                  maxquant_row_ids = paste0(evidence_id, collapse=";")) %>%
        ungroup() %>%
        mutate(replicate = toupper(replicate))
      output_tbl <- temp_tbl
      if(!is.na(additional_cols) & additional_cols != "") {
        output_tbl <- temp_tbl %>%
          mutate_at(additional_cols, toupper)
      }
      return(output_tbl)
    }
    summary_funcs <- list(mean=mean, median=median, sum=sum)
    summarised_long_tbl_list <- purrr::map(summary_funcs, 
                                           ~group_summary(all_phos_sites_long_tbl, additional_cols, .))
    return( summarised_long_tbl_list)
  }
  
  summarised_long_tbl_list <- uniquePhosphositesSummariseLongList(paralog_sites_long,
                                                                  additional_cols)
  
  #----- Summarise the abundance values for each unique phosphosites (mean, median, sum), return table in wide format
  uniquePhosphositesSummariseWideList <- function(summarised_long_tbl_list,
                                                  additional_cols=c("experiment")) {
    cols_to_use <- c("replicate")
    summarised_wide_tbl_list_edited <- NA
    
    if(!is.na( additional_cols) & additional_cols != "") {
      cols_to_use <- c("replicate", additional_cols)
      experiment_col <- additional_cols[[1]]
      summarised_wide_tbl_list_edited <- 
        purrr::map(summarised_long_tbl_list, 
                   function(input_table){
                     output_table <- input_table %>% 
                       mutate(maxquant_row_ids = paste0(paste(!!rlang::sym(experiment_col), 
                                                              sep="_"), 
                                                        "(",
                                                        maxquant_row_ids, 
                                                        ")")) %>% 
                       pivot_wider(id_cols = c("uniprot_acc", "gene_names", "protein_site_positions", 
                                               "phos_15mer_seq", "maxquant_row_ids"), 
                                   names_from = all_of(cols_to_use), 
                                   values_from=c("value")) %>% 
                       unite("sites_id", uniprot_acc, gene_names, protein_site_positions, 
                             phos_15mer_seq, sep="!") 
                     return(output_table)})
    } else {
      summarised_wide_tbl_list_edited <- 
        purrr::map(summarised_long_tbl_list, 
                   function(input_table){ 
                     output_table <- input_table %>% 
                       pivot_wider(id_cols = c("uniprot_acc", "gene_names", "protein_site_positions", 
                                               "phos_15mer_seq", "maxquant_row_ids"), 
                                   names_from = all_of(cols_to_use), 
                                   values_from=c("value")) %>% 
                       unite("sites_id", uniprot_acc, gene_names, protein_site_positions, 
                             phos_15mer_seq, sep="!") 
                     return(output_table)})
    }
    
    #- Summarize MaxQuant evidence IDs from different multiplex experiment
    clean_maxquant_ids <- function(input_tab ) {
      maxquant_ids_tbl <- input_tab %>%
        group_by(sites_id) %>%
        summarise(maxquant_row_ids = paste(maxquant_row_ids, collapse=",")) %>%
        ungroup()
      values_tbl <- input_tab %>%
        dplyr::select(-maxquant_row_ids) %>%
        group_by(sites_id) %>%
        summarise_all(~sum(., na.rm=TRUE)) %>%
        ungroup()
      output_tab <- values_tbl %>%
        left_join(maxquant_ids_tbl, by="sites_id") %>%
        relocate(maxquant_row_ids, .after="sites_id")
      colnames(output_tab) <- colnames(output_tab) %>%
        toupper() %>%
        str_replace_all("SITES_ID", "sites_id")  %>%
        str_replace_all("MAXQUANT_ROW_IDS", "maxquant_row_ids") 
      return( output_tab)
    }
    summarised_wide_tbl_cln_list <- purrr::map( summarised_wide_tbl_list_edited, clean_maxquant_ids)
    return( summarised_wide_tbl_cln_list)
  }
  
  summarised_wide_tbl_list <- uniquePhosphositesSummariseWideList(summarised_long_tbl_list,
                                                                  additional_cols)
  
  #----- Output files
  results_list <- list(summarised_wide_list = summarised_wide_tbl_list, 
                       summarised_long_list = summarised_long_tbl_list,
                       all_phos_sites_wide  = all_phos_sites_wide_tbl,
                       all_phos_sites_long  = all_phos_sites_long_tbl )
  
  output_files <- map_chr(list("mean", "median", "sum"),
                          ~file.path(output_dir, paste0(., "_phoshpsites.tsv")))
  
  walk2(results_list$summarised_wide_list,
        output_files,
        ~vroom::vroom_write(.x, .y))
  
  vroom::vroom_write(results_list$all_phos_sites_wide,
                     file.path(output_dir, 
                               "all_phos_sites_wide_tbl.tsv"))
  
  saveRDS(results_list$summarised_long_list,
          file.path(output_dir, 
                    "summarised_long_tbl_list.RDS"))
  
  saveRDS(results_list$summarised_wide_list,
          file.path(output_dir, 
                    "summarised_wide_tbl_list.RDS"))
  
  sample_names_file <- file.path(output_dir, "sample_names.tab")
  sample_names <- colnames(results_list$summarised_wide_list[[1]])[c(-1,-2)]
  vroom::vroom_write(data.frame(sample_names = t(t(sample_names))), 
                     sample_names_file)
  
  #- Fix Column names
  # colnames(evidence_tbl_filt) <- gsub(pattern = "REPORTER_INTENSITY_CORRECTED_", 
  #                                     replacement = "", 
  #                                     x = colnames(evidence_tbl_filt))
}