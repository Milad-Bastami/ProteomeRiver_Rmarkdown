#---------- Phosphoproteome Kinswingr
#----- Packages
# suppressPackageStartupMessages({
#   p_load(UniProt.ws)
#   p_load(BiocParallel)
#   p_load(KinSwingR)
#   p_load(RColorBrewer)
#   p_load(future)
# })

#----- Variables
# #- Working Directory
# setwd("/home/ubuntu/BMP_Block_LM/20240213_Proteomeriver_Rmd_BIO_25_LM")
# #- Directory path for temporary files
# tmp_dir <-"results/phosphoproteomics/cache"
# #- The number of cores used for the computation
# num_cores <- 1
# #- Random seed
# random_seed <- 123456
# #- The number iterations for scoring each substrate against kinase motifs
# motif_score_iteration <- 1000
# #- The number iterations for swing score calcualtions
# swing_iteration <- 1000
# #- The minimum number of known substrates for each kinase
# min_num_sites_per_kinase <- 10
# #- A file that contains the dictionary to convert uniprot accession to gene symbol. 
# #- Uses the column specified in 'protein_id_lookup_column' flag for protein ID. 
# #- Uses the column specified in 'gene_symbol_column' for the gene symbol column
# uniprot_to_gene_symbol_file <- "data/HomoSapiens/uniprot/data.tab"
# #- The name of the column that contained the protein ID to convert into gene symobl
# protein_id_lookup_column <- "Entry"
# #- The name of the column that contained the gene symobl
# gene_symbol_column <- "Gene Names"
# #- The p-value cutoff for identifying a kinase as significant by KinSwingR
# p_value_cutoff <- 0.05
# #- The p-value cutoff for showing the kinase name in the volcano plot
# ggrepel_p_value_cutoff <- 0.1
# #- The directory in which PhosphositePlus data is stored
# phosphosite_db_dir <- "data/HomoSapiens/phosphosites"
# #- The path to the UniProt tab-separated data file
# uniprot_kinase_file <- "data/HomoSapiens/uniprot/pkinfam.tab"
# #- Results table in which the phosphorylation log fold-change is normalized 
# #- by protein log fold-change
# norm_phos_logfc_file <- "results/phosphoproteomics/norm_phos_by_prot_abundance/norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv" 
# #- File listing the UniProt accession of atypical and other kinases and their UniProt Keywords
# uniprot_other_kinase_file <- "atypical_and_other_kinases_data.RDS"
# #- Column name in the input file that contains the log fold-change values of the phosphosites
# log_fc_column_name <- "norm_phos_logFC"
# #- The NCBI taxonomy ID of the organism being investigated (M.musculus=10090, H.sapien=9606)
# taxonomy_id <- 9606
# #- Column name in the input file that contains the false discovery rate values of the phosphosites
# fdr_column_name <- "combined_q_mod"
# #- If NA, include both single-site and multi-site phosphorylation. If TRUE, include 
# #- only multi-site phosphorylation. If FALSE, include only single-site phosphorylation
# is_multisite <- NA
# #- If TRUE, reuse result from files saved from previous run
# reuse_old <- FALSE
# #- A comma separated strings to indicate the fortmat output for the plots [pdf,png,svg]
# plots_format <- list("pdf","png")

#----- Phosphoproteome Kinswingr ST
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phos_kinswingr_ST"
# #- Specificity of the kinase. One of Ser/Thr kinase = ST, Tyr kinase = Y, or 
# #- Ser/Thr/Tyr kinase = STY
# kinase_specificity <- "ST"

#----- Phosphoproteome Kinswingr STY
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phos_kinswingr_STY"
# #- Specificity of the kinase. One of Ser/Thr kinase = ST, Tyr kinase = Y, or 
# #- Ser/Thr/Tyr kinase = STY
# kinase_specificity <- "STY"

#----- Phosphoproteome Kinswingr STY
# #- Directory path for all results files
# output_dir <- "results/phosphoproteomics/phos_kinswingr_Y"
# #- Specificity of the kinase. One of Ser/Thr kinase = ST, Tyr kinase = Y, or 
# #- Ser/Thr/Tyr kinase = STY
# kinase_specificity <- "Y"

#----- Function
getPhosphoproteomeKinswingr <- function(
    tmp_dir,
    num_cores,
    random_seed,
    motif_score_iteration,
    swing_iteration,
    min_num_sites_per_kinase,
    uniprot_to_gene_symbol_file,
    protein_id_lookup_column,
    gene_symbol_column,
    p_value_cutoff,
    ggrepel_p_value_cutoff,
    phosphosite_db_dir,
    uniprot_kinase_file,
    norm_phos_logfc_file,
    uniprot_other_kinase_file,
    log_fc_column_name,
    taxonomy_id,
    fdr_column_name,
    is_multisite,
    reuse_old,
    plots_format,
    output_dir,
    kinase_specificity
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
  
  #----- Reading UniProt kinases list file
  uniprot_kinase_tbl <- vroom::vroom(uniprot_kinase_file)
  
  #----- Reading UniProt data file
  if(file.exists(uniprot_to_gene_symbol_file)) {
    uniprot_tab_delimited_tbl <- vroom::vroom(file.path(uniprot_to_gene_symbol_file))
  }
  
  #----- Download information from UniProt
  #- Clean Isoform Number
  #- clean_isoform_number("Q8K4R4-2")
  cleanIsoformNumber <- function(string) {
    str_replace(string, "-\\d+$", "")
  }
  
  #- The UniProt.ws::select function limits the number of keys queried to 100. 
  #- This gives a batch number for it to be queried in batches
  batchQueryEvidenceHelper <- function(uniprot_acc_tbl, 
                                       uniprot_acc_column){
    #- 100 is the maximum number of queries at one time
    all_uniprot_acc <- uniprot_acc_tbl |>
      dplyr::select({{uniprot_acc_column}}) |>
      mutate(Proteins = str_split({{uniprot_acc_column}}, ";")) |>
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
  
  uniprot_other_kinase_file <- file.path(tmp_dir, 
                                         uniprot_other_kinase_file)
  
  if(!file.exists(uniprot_other_kinase_file)) {
    human_taxonomy_id <- 9606
    up <- UniProt.ws(taxId=human_taxonomy_id)
    list_of_sp_columns <- c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", 
                            "PROTEIN-NAMES", "LENGTH", "ENSEMBL", "GO-ID", 
                            "KEYWORDS", "protein_existence", "annotation_score", 
                            "reviewed", "gene_names", "protein_name", "length", 
                            "xref_ensembl", "go_id", "keyword")
    up_cls <- unlist(columns(up))
    list_intersect <- intersect(list_of_sp_columns, up_cls)
    
    # if(length(setdiff(list_of_sp_columns, list_intersect)) > 0) {
    #   print("UniProt fields not found: %s", 
    #         paste(setdiff(list_of_sp_columns, list_intersect), sep=", "))
    # }
    
    my_keytype <- "UniProtKB"
    if("UNIPROTKB" %in% keytypes(up)) {
      my_keytype <- "UNIPROTKB"
    }
    
    uniprot_acc_tbl <- uniprot_kinase_tbl |>
      dplyr::filter(Family == "Other" | 
                      str_detect(Family, "Atypical")) |>
      dplyr::distinct(uniprot_acc_human) |>
      dplyr::rename(uniprot_acc = uniprot_acc_human)
    
    atypical_and_other <- batchQueryEvidence(uniprot_acc_tbl, 
                                             uniprot_acc_column = "uniprot_acc", 
                                             uniprot_handle = up,
                                             uniprot_columns = list_intersect, 
                                             uniprot_keytype = my_keytype)
    
    saveRDS(atypical_and_other, 
            file.path(uniprot_other_kinase_file))
  }
  
  atypical_and_other <- readRDS(uniprot_other_kinase_file)
  
  #----- Read PhosphoSitePlus (PSP) kinase-substrate table
  ks_file <- file.path(phosphosite_db_dir, 
                       "Kinase_Substrate_Dataset")
  
  ks_tbl <- vroom::vroom(ks_file, skip=3) |>
    mutate(SUB_MOD_RSD_CLN = str_replace_all(SUB_MOD_RSD, 
                                             "([A-Z])(\\d+)", 
                                             "\\1 \\2")) |>
    separate(SUB_MOD_RSD_CLN, 
             into=c("residue", "position"))
  
  #----- Set number of cores for parallel computation
  register(SnowParam(workers = num_cores))
  
  #----- Check log fold-change and FDR column names exists in input file
  # if(!log_fc_column_name %in% colnames(de_phos)) {
  #   print(paste0("Column '%s' is not found in the input table.", 
  #                log_fc_column_name))
  # }
  # if(fdr_column_name %in% colnames(de_phos)) {
  #   print(paste0("Column '%s' is not found in the input table.", 
  #                fdr_column_name))
  # }
  # if(!log_fc_column_name %in% colnames(de_phos) | 
  #    !fdr_column_name %in% colnames(de_phos)) {
  #   stop()
  # }
  
  #----- Filter the correct subset of kinases for the analysis
  uniprot_kinases <- NA
  
  if("UNIPROTKB" %in% colnames(atypical_and_other)) {
    atypical_and_other <- atypical_and_other |>
      rename(Keywords = "KEYWORDS", 
             Entry = "UNIPROTKB")
  }
  
  if(kinase_specificity == "ST") {
    uniprot_kinases <- uniprot_kinase_tbl |>
      dplyr::filter(str_detect(Family, "Ser/Thr") | 
                      str_detect(Family, "Atypical") | 
                      str_detect(Family, "Other")) |> 
      left_join(atypical_and_other |> 
                  dplyr::select( one_of( "Entry", "Keywords")),
                by = c("uniprot_acc_human" = "Entry")) |>
      dplyr::filter(str_detect(Family, "Ser/Thr") | 
                      (str_detect(Keywords, "Serine/threonine-protein kinase") & 
                         (!str_detect(Keywords, "Tyrosine-protein kinase"))))
  } else if(kinase_specificity == "Y") {
    uniprot_kinases <- uniprot_kinase_tbl |>
      dplyr::filter(str_detect(Family, "Tyr") | 
                      str_detect(Family, "Atypical") | 
                      str_detect(Family, "Other")) |>
      left_join(atypical_and_other |>
                  dplyr::select(one_of( "Entry", "Keywords")),
                by = c("uniprot_acc_human" = "Entry")) |>
      dplyr::filter(str_detect(Family, "Tyr") |
                      ((!str_detect(Keywords, "Serine/threonine-protein kinase")) & 
                         str_detect(Keywords, "Tyrosine-protein kinase")))
  } else if(kinase_specificity == "STY") {
    uniprot_kinases <- uniprot_kinase_tbl |>
      dplyr::filter(str_detect(Family, "Atypical") | 
                      str_detect(Family, "Other")) |>
      left_join(atypical_and_other |>
                  dplyr::select(one_of("Entry", "Keywords")),
                by = c("uniprot_acc_human" = "Entry")) |>
      dplyr::filter(str_detect(Keywords, "Serine/threonine-protein kinase") & 
                      str_detect(Keywords, "Tyrosine-protein kinase"))
  }
  
  #----- Filter the correct subset of substrates for the analysis
  phosphositeplus <- ks_tbl |>
    dplyr::select(GENE, KINASE, `SITE_+/-7_AA`, SUB_MOD_RSD) |>
    distinct() |>
    dplyr::rename(kinase_class = "KINASE",
                  substrate = "SITE_+/-7_AA",
                  kinase_gene = "GENE") |>
    dplyr::mutate(substrate = toupper(substrate)) |>
    dplyr::mutate(kinase_class = str_replace(kinase_class, " ", "_")) |>
    arrange(kinase_gene, kinase_class, substrate) |>
    mutate(kinase_gene = toupper(kinase_gene),
           kinase_class = toupper(kinase_class),
           residue = purrr::map_chr(SUB_MOD_RSD, 
                                    \(x){str_sub(x, 1,1 )})) |>
    dplyr::select(-SUB_MOD_RSD)
  
  kinases_gene_to_substrate_counts <- phosphositeplus |>
    left_join(uniprot_kinases |>
                dplyr::select(gene_name) |>
                dplyr::mutate(is_uniprot_a = 1),
              by = c("kinase_gene" = "gene_name")) |>
    left_join(uniprot_kinases |>
                dplyr::select(gene_name) |>
                dplyr::mutate(is_uniprot_b = 1),
              by = c("kinase_class" = "gene_name")) |>
    dplyr::filter(is_uniprot_a == 1 | is_uniprot_b == 1) |>
    dplyr::select(-is_uniprot_a, -is_uniprot_b) |>
    dplyr::filter(str_detect(kinase_specificity, residue)) |>
    group_by(kinase_gene) |>
    summarise(counts = n()) |>
    ungroup() 
  
  kinases_to_include <- kinases_gene_to_substrate_counts |>
    dplyr::filter(counts >= min_num_sites_per_kinase)
  
  phosphositeplus_filt <- phosphositeplus |>
    inner_join(kinases_to_include, 
               by = "kinase_gene") |>
    dplyr::filter(str_detect(kinase_specificity, residue)) |>
    dplyr::select(-counts, -kinase_class, -residue) |>
    distinct() |>
    as.matrix()
  
  #----- Build the position-specific scoring matrices (PSSM)
  pwms <- buildPWM(as.matrix(phosphositeplus_filt))
  
  #----- Clean up phoshopsite log fold-change table
  annotated_data_pre_residue_filter <- de_phos |>
    mutate(peptide = sequence) |>
    mutate(positions_peptide = position) |>
    mutate(uniprot_acc = str_split(uniprot_acc, ":") |> 
             purrr::map_chr(1)) |>
    mutate(gene_name = str_split(gene_name, ":") |> 
             purrr::map_chr(1))  |>
    mutate(position = str_split(position, ":") |> 
             purrr::map_chr(1)) |>
    mutate(peptide =  str_split(peptide, ":") |> 
             purrr::map_chr(1)) |>
    mutate(residue =  str_split(residue, ":") |> 
             purrr::map_chr(1)) |>
    mutate(position = str_split(position, "\\|") |> 
             purrr::map_chr(1) |> 
             str_replace_all("\\(|\\)", "")) |>
    mutate(is_multisite = case_when(str_detect(position, ";") ~ TRUE,
                                    TRUE ~ FALSE)) |>
    separate_rows(peptide, position, residue, sep=";") |>
    dplyr::mutate(peptide_copy = peptide) |>
    dplyr::filter(!str_detect(peptide, "X")) |>
    dplyr::select(sites_id, comparison, uniprot_acc, gene_name, 
                  peptide, peptide_copy, position, residue,
                  one_of(c(as.character(log_fc_column_name), as.character(fdr_column_name))), 
                  is_multisite) |>
    dplyr::filter(str_detect(kinase_specificity, residue)) |>
    dplyr::filter(!is.na(peptide)) |>
    dplyr::filter(str_sub( peptide, 8, 8) == residue) |>
    dplyr::rename(fc = log_fc_column_name,
                  pval = fdr_column_name)
  
  #----- Filtering single-site or multisite phosphorylation
  if(!is.na(is_multisite)) {
    annotated_data_multisite_filtered <- annotated_data_pre_residue_filter |>
      dplyr::filter(is_multisite == is_multisite)
    if(is_multisite) {
      print("Keeping only multi-site phosphorylation.")
    } else {
      print("Keeping only single-site phosphorylation.")
    }
  } else {
    print("Keeping both single-site and multi-site phosphorylation.")
    annotated_data_multisite_filtered <- annotated_data_pre_residue_filter
  }
  
  #----- For sites with multisites information, get best p-value (and log fc if there are ties) 
  #- for each comparison, UniProt accession, and position combination")
  best_p_value <- annotated_data_multisite_filtered |>
    group_by(uniprot_acc, position, comparison) |>
    dplyr::summarise(best_p_val = min(pval)) |>
    ungroup()
  
  best_log_fc_pval_site_join <- c("uniprot_acc", "position", "best_p_val", "comparison")
  names( best_log_fc_pval_site_join) <- c("uniprot_acc", "position", "pval", "comparison")
  
  best_logfc_pval_site <- annotated_data_pre_residue_filter |>
    inner_join(best_p_value, 
               by = best_log_fc_pval_site_join) |>
    group_by(uniprot_acc, position, comparison, pval) |>
    dplyr::summarise(best_abs_log_fc =  max(abs(fc))) |>
    ungroup()
  
  best_site <- annotated_data_pre_residue_filter |>
    mutate(abs_log_fc = abs(fc)) |>
    inner_join(best_logfc_pval_site, 
               by = c("uniprot_acc" = "uniprot_acc",
                      "position" = "position",
                      "pval" = "pval",
                      "abs_log_fc" = "best_abs_log_fc",
                      "comparison" = "comparison") ) |>
    dplyr::select(-abs_log_fc)
  
  annotated_data <- best_site |>
    unite(annotation, uniprot_acc, gene_name, position, peptide_copy , sep="|") |>
    distinct()
  
  #----- Perform analysis of data from each contrast
  grouped_annotated_data <- annotated_data |>
    dplyr::select(-is_multisite, -sites_id, -residue) |>
    distinct() |>
    group_by(comparison) |>
    nest() |>
    ungroup()
  
  #----- Score each phosphorylation site for similarity to kinase-specific motif
  rds_dir <- output_dir
  # if (no_backup) {
  #   rds_dir <- paste(output_dir, "_prev", sep = "")
  # }
  
  if(reuse_old == TRUE & 
     file.exists(file.path(rds_dir, 
                           paste0("kinswingr_scores_list_", 
                                  kinase_specificity, 
                                  ".RDS")))) {
    
    scores_list <- readRDS(file.path(rds_dir, 
                                     paste0("kinswingr_scores_list_", 
                                            kinase_specificity, 
                                            ".RDS")))
    
    copy_attempt <- file.copy(from = file.path(rds_dir, 
                                               paste0("kinswingr_scores_list_", 
                                                      kinase_specificity, 
                                                      ".RDS")),
                              to =  file.path(output_dir, 
                                              paste0("kinswingr_scores_list_", 
                                                     kinase_specificity, 
                                                     ".RDS")))
  } else {
    # set seed for reproducible results
    set.seed(random_seed)
    
    print(paste("Number of sites per contrast: ",
                paste(names(grouped_annotated_data$data), 
                      collapse=", "),
                paste(purrr::map_int(grouped_annotated_data$data, 
                                     \(x){nrow(as.data.frame(x))}), 
                      collapse=", ")))
    
    ## Use the number of data to determine the number of iterations
    num_of_iterations <- purrr::map_dbl(grouped_annotated_data$data,
                                        \(x){ 
                                          num_of_rows <- nrow(as.data.frame(x))
                                          if(num_of_rows < motif_score_iteration) {
                                            motif_score_iteration_to_use <- num_of_rows - (num_of_rows %% 10)
                                          } else {
                                            motif_score_iteration
                                          }
                                        })
    
    scores_list <-  grouped_annotated_data |>
      dplyr::mutate(pwms_scores = purrr::map2(data, 
                                              num_of_iterations, 
                                              \(.x, .y){
                                                scoreSequences(input_data = as.data.frame(.x),
                                                               pwm_in = pwms,
                                                               n = .y )
                                              }))
    
    saveRDS(scores_list, 
            file.path(output_dir, paste0("kinswingr_scores_list_", kinase_specificity, ".RDS")))
  }
  
  purrr::walk2(scores_list$data,
               scores_list$comparison,
               \(.x, .y){
                 vroom::vroom_write(.x |> 
                                      mutate(comparison = .y), 
                                    file.path(output_dir,
                                              paste0("input_data_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv")))
               })
  purrr::walk2(scores_list$pwms_scores,
               scores_list$comparison,
               \(.x, .y){
                 vroom::vroom_write(.x$peptide_scores |>
                                      mutate( comparison = .y),
                                    file.path(output_dir,
                                              paste0("peptide_scores_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv")))
               })
  purrr::walk2(scores_list$pwms_scores,
               scores_list$comparison,
               \(.x, .y){ 
                 vroom::vroom_write(.x$peptide_p |>
                                      mutate(comparison = .y),
                                    file.path(output_dir,
                                              paste0("peptide_p_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv"))) 
               })
  purrr::walk2(scores_list$pwms_scores,
               scores_list$comparison,
               \(.x, .y){ 
                 vroom::vroom_write(.x$background |>
                                      mutate(comparison = .y),
                                    file.path(output_dir,
                                              paste0("peptide_background_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv")))
               })
  
  #----- Perform randomization analysis with Swing
  # set seed for reproducible results
  set.seed(random_seed)
  
  if(reuse_old == TRUE & 
     file.exists(file.path(rds_dir, 
                           paste0("kinswingr_swing_out_list_", 
                                  kinase_specificity, 
                                  ".RDS")))) {
    
    swing_out_list <- readRDS(file.path(rds_dir,
                                        paste0("kinswingr_swing_out_list_", 
                                               kinase_specificity, 
                                               ".RDS")))
    
    copy_attempt <- file.copy(from = file.path(rds_dir, 
                                               paste0("kinswingr_swing_out_list_", 
                                                      kinase_specificity, 
                                                      ".RDS")),
                              to = file.path(output_dir, 
                                             paste0("kinswingr_swing_out_list_", 
                                                    kinase_specificity, 
                                                    ".RDS")))
  } else {
    swing_out_list <- scores_list |>
      dplyr::mutate(swing_result = purrr::map2(data, 
                                               pwms_scores, 
                                               \(.x, .y) swing(input_data = as.data.frame(.x),
                                                               pwm_in = pwms,
                                                               pwm_scores = .y,
                                                               permutations = swing_iteration,
                                                               return_network = TRUE)))
    #- This will produce two tables, one is a network for use with e.g. Cytoscape and 
    #- the other is the scores. To access the scores:
    head(swing_out_list$swing_result[[1]]$scores)
    saveRDS(swing_out_list, file.path(output_dir,
                                      paste0("kinswingr_swing_out_list_", 
                                             kinase_specificity, ".RDS")))
  }
  # swing_out_list<- readRDS(file.path(output_dir,
  #                                    paste0("kinswingr_swing_out_list_", 
  #                                           kinase_specificity, ".RDS")))
  
  purrr::walk2(swing_out_list$swing_result,
               swing_out_list$comparison,
               \(.x, .y){ 
                 vroom::vroom_write(.x$scores |>
                                      mutate(comparison = .y),
                                    file.path(output_dir,
                                              paste0("KinSwingR_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv")))
               })
  purrr::walk2(swing_out_list$swing_result,
               swing_out_list$comparison,
               \(.x, .y){ 
                 vroom::vroom_write(.x$network |>
                                      mutate(comparison = .y),
                                    file.path(output_dir,
                                              paste0("network_", 
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv"))) 
               })
  
  #----- Plot kinase-level volcano plot
  kinase_volcano_plot <- function(
    x_label_name,
    kinase_type = "ST",
    list_position,
    brewer_pal_set = "Set1",
    p_value_cutoff = 0.05,
    ggrepel_p_value_cutoff = 0.05,
    n_cutoff = 10) {
    
    swing_scores_tbl <- swing_out_list$swing_result[[list_position]]$scores
    
    kinase_type_string <- case_when(kinase_type == "ST" ~ "predicted Ser/Thr PK activity",
                                    kinase_type == "Y" ~ "predicted Tyr PK activity",
                                    kinase_type == "STY" ~ "predicted Ser/Thr/Tyr PK activity")
    
    my_colours <- setdiff(brewer.pal(length(swing_out_list$comparison)+1, 
                                     brewer_pal_set), 
                          "#FFFF33")
    
    significance_tbl <- swing_scores_tbl |>
      mutate(p_value = case_when(p_greater < p_less ~ p_greater,
                                 p_less <= p_greater ~ p_less)) |>
      mutate(colour = case_when(p_value < p_value_cutoff ~ "Significant",
                                TRUE ~ "Not significant")) |>
      mutate(colour = factor(colour, 
                             levels=c("Not significant", "Significant"))) |>
      mutate(colour_vector = case_when(colour == "Not significant" ~ "grey",
                                       colour == "Significant" ~ my_colours[list_position]))
    
    scale_colour_vector <- significance_tbl |>
      arrange(colour) |>
      distinct(colour_vector) |>
      pull(colour_vector)
    
    volcano_plot_dat <- significance_tbl |>
      mutate(kinase_label = case_when(p_value < ggrepel_p_value_cutoff & 
                                        n > n_cutoff ~ kinase,
                                      TRUE ~ ""))
    
    if(nrow(volcano_plot_dat) > 0 ) { 
      volcano_plot <- volcano_plot_dat |>
        ggplot(aes(x = swing, y = -log10(p_value + .Machine$double.eps), 
                   size = n, col = colour, label = kinase_label, alpha = 0.5)) +
        theme_bw() +
        geom_point() +
        scale_color_manual(values = scale_colour_vector, 
                           name = "Is significant" ) +
        scale_size_continuous(name = "Num. substrates") +
        ylab(expression("Significance, -log"[10]*"(P)")) +
        xlab(paste0(x_label_name[list_position], ", ", kinase_type_string)) +
        geom_hline(aes(yintercept = -log10(p_value_cutoff)), linetype=3) +
        geom_text_repel(show.legend = FALSE, size = 5) +
        theme(text = element_text(size = 20))
      return(volcano_plot)
    } else {
      return(NA)
    }
  }
  
  partial_kinase_volcano_plot <- partial(kinase_volcano_plot, 
                                         x_label_name = swing_out_list$comparison,
                                         kinase_type = kinase_specificity,
                                         p_value_cutoff = p_value_cutoff,
                                         n_cutoff = min_num_sites_per_kinase,
                                         ggrepel_p_value_cutoff = ggrepel_p_value_cutoff)
  
  kinase_volcano_plot <- purrr::map(seq_along(swing_out_list$comparison),
                                    \(x){partial_kinase_volcano_plot(list_position = x)})
  
  for(format_ext in plots_format) {
    purrr::walk2(kinase_volcano_plot, swing_out_list$comparison,
                 \(.x, .y) { 
                   if(length(.x) > 1 ) { 
                     ggsave(filename = file.path(output_dir,
                                                 paste0("kinase_volcano_plot_",
                                                        kinase_specificity, 
                                                        "_", .y, ".", format_ext)), 
                            plot=.x)}})
  }
  
  #----- Plot Average logFC Known Sites
  plotAvgLogfcKnownSites <- function(list_position,
                                     y_label_name,
                                     brewer_pal_set= "Set1",
                                     kinase_type,
                                     p_value_cutoff = 0.05,
                                     ggrepel_p_value_cutoff = 0.05) {
    
    kinase_type_string <- case_when(kinase_type == "ST" ~ "predicted Ser/Thr PK activity",
                                    kinase_type == "Y" ~ "predicted Tyr PK activity",
                                    kinase_type == "STY" ~ "predicted Ser/Thr/Tyr PK activity")
    
    my_colours <- brewer.pal(length(swing_out_list$comparison), brewer_pal_set)
    
    my_comparison <- swing_out_list$comparison[[list_position]]
    
    fc_pval_tab <- scores_list$data[[list_position]]
    
    up_or_down_kinases <- swing_out_list$swing_result[[list_position]]$scores |>
      dplyr::filter(p_greater < p_value_cutoff | p_less < p_value_cutoff)
    
    motif_score_table <- scores_list$pwms_scores[[list_position]]$peptide_scores |>
      pivot_longer(cols = !(contains("annotation") | contains("peptide")),
                   names_to = "kinase",
                   values_to = "motif.score")
    
    peptide_p_long_tbl <- scores_list$pwms_scores[[list_position]]$peptide_p |>
      pivot_longer(cols = !(contains("annotation") | contains("peptide")),
                   names_to = "kinase",
                   values_to = "motif.p.value")
    
    selected_scores_list <- peptide_p_long_tbl |>
      left_join(motif_score_table, 
                by = c("annotation" = "annotation",
                       "peptide" = "peptide",
                       "kinase" = "kinase")) |>
      inner_join(up_or_down_kinases, 
                 by = c( "kinase" = "kinase")) |>
      inner_join(fc_pval_tab, 
                 by = c("annotation" = "annotation",
                        "peptide" = "peptide")) |>
      left_join(annotated_data |>
                  dplyr::select(annotation, sites_id),
                by = c("annotation" = "annotation")) |>
      left_join(de_phos |>
                  dplyr::select(sites_id, KINASE), 
                by = c("sites_id" = "sites_id")) |>
      distinct() |>
      separate_rows(KINASE, sep= "//") |>
      dplyr::mutate(KINASE = toupper(KINASE) ) |>
      dplyr::filter(kinase == KINASE)
    
    avg_logfc_vs_swing_score <- selected_scores_list |>
      group_by(kinase, swing, p_greater, p_less) |>
      summarise(avg_fc = mean(fc),
                median_fc = median(fc),
                total_fc = sum(fc),
                counts = n()) |>
      ungroup() |>
      dplyr::filter(counts >= 5 )
    
    avg_logfc_vs_swing_tbl <- avg_logfc_vs_swing_score |>
      dplyr::mutate(p_value = case_when(swing > 0 ~ p_greater,
                                        swing < 0 ~ p_less,
                                        TRUE ~ 1))  |>
      mutate(colour = case_when((p_value < p_value_cutoff |
                                   p_value < p_value_cutoff ) ~ "Significant",
                                TRUE ~ "Not significant")) |>
      dplyr::mutate(my_kinase = case_when(counts >= 10 |
                                            (p_greater < ggrepel_p_value_cutoff | 
                                               p_less < ggrepel_p_value_cutoff) ~ kinase,
                                          TRUE ~ ""))
    
    if(nrow(avg_logfc_vs_swing_tbl) > 0) {
      avg_logfc_vs_swing_score_plot <- avg_logfc_vs_swing_tbl |>
        ggplot(aes(y = avg_fc, x = swing, label = my_kinase, 
                   size = counts, color = colour, alpha = 0.5 )) +
        geom_point() +
        geom_text_repel(show.legend = FALSE, size = 5) +
        scale_color_manual(values = c("grey", my_colours[list_position]), name="Is significant") +
        ylab(substitute(paste(a, ", Mean substrate log"[2], "(intensity) difference"), 
                        list(a=y_label_name[list_position]))) +
        xlab(paste0(y_label_name[list_position], ", ", kinase_type_string)) +
        theme(text = element_text(size = 20))
      return(avg_logfc_vs_swing_score_plot)
    }
    return(NA)
  }
  
  # scores_list <- readRDS( "/home/ignatius/PostDoc/2021/Tautomycetin_2021/Results/phos_kinswingr_ST_exp1/kinswingr_scores_list_ST.RDS" )
  # swing_out_list <- readRDS( "/home/ignatius/PostDoc/2021/Tautomycetin_2021/Results/phos_kinswingr_ST_exp1/kinswingr_swing_out_list_ST.RDS" )
  
  ## Do this only if the species is Human or Mouse
  if(taxonomy_id %in% c(9606, 10090)) {
    print("Plot average log fold-change of known phoshosites.")
    
    partial_plotAvgLogfcKnownSites <- partial(plotAvgLogfcKnownSites,
                                              y_label_name = swing_out_list$comparison,
                                              brewer_pal_set= "Set1",
                                              kinase_type = kinase_specificity,
                                              p_value_cutoff = p_value_cutoff,
                                              ggrepel_p_value_cutoff = ggrepel_p_value_cutoff)
    
    kinase_avg_logfc <- purrr::map(seq_along(swing_out_list$comparison),
                                   \(x) partial_plotAvgLogfcKnownSites(x) )
    
    for(format_ext in plots_format) {
      purrr::walk2(kinase_avg_logfc,
                   swing_out_list$comparison,
                   \(.x, .y){ if ( length(.x) > 1) { 
                     ggsave(filename = file.path(output_dir,
                                                 paste0("kinase_avg_logfc_of_known_sites_",
                                                        kinase_specificity, 
                                                        "_", .y, ".", format_ext )),
                            plot=.x) }})
    }
  }
  
  #----- Compile KinSwingR results 
  list_of_kinswinger_columns <- c("substrate_uniprot_acc",
                                  "subsrate_gene_symbol",
                                  "phosphosite_position",
                                  "sequence_context",
                                  "kinase_uniprot_acc",
                                  "kinase_gene_symbol",
                                  "kinsae_gene_name_uniprot",
                                  "kinase_family",
                                  "phosphosite_log2FC",
                                  "phosphosite_fdr_value",
                                  "motif.score",
                                  "motif.p.value",
                                  "swing_kinswingr",
                                  "p_value_kinswingr",
                                  "pos_kinswingr",
                                  "neg_kinswingr",
                                  "all_kinswingr",
                                  "pk_kinswingr",
                                  "nk_kinswingr",
                                  "swing_raw_kinswingr",
                                  "n_kinswingr",
                                  "p_greater_kinswingr",
                                  "p_less_kinswingr",
                                  "is_kinase_phosphorylated",
                                  "known_upstream_kinase",
                                  "prediction_match_known_kinase",
                                  "reactome_term",
                                  "ON_FUNCTION",
                                  "ON_PROCESS",
                                  "ON_PROT_INTERACT",
                                  "ON_OTHER_INTERACT",
                                  "REG_SITES_NOTES",
                                  "kinase_uniprot_id_human",
                                  "kinase_uniprot_acc_human",
                                  "kinase_uniprot_id_mouse",
                                  "kinase_uniprot_acc_mouse",
                                  "sites_id",
                                  "substrate_name")
  
  compileKinswingerResults <- function(list_position) {
    
    my_comparison <- swing_out_list$comparison[[list_position]]
    
    fc_pval_tab <- scores_list$data[[list_position]]
    
    up_or_down_kinases <- swing_out_list$swing_result[[list_position]]$scores |>
      dplyr::filter(p_greater < 0.2 | p_less < 0.2) |>
      dplyr::rename(kinase_gene ="kinase")
    
    motif_score_table <- scores_list$pwms_scores[[list_position]]$peptide_scores |>
      pivot_longer(cols = !(contains("annotation") | contains("peptide")),
                   names_to ="kinase_gene",
                   values_to = "motif.score")
    
    is_phosphoylated_tbl <- annotated_data |>
      dplyr::filter(pval < p_value_cutoff) |>
      dplyr::mutate(uniprot_acc = str_split(annotation, "\\|") |> 
                      purrr::map_chr(1)) |>
      distinct(uniprot_acc) |>
      mutate(is_kinase_phosphorylated = 1)
    
    selected_columns <- intersect(colnames(de_phos),
                                  c("sites_id", "PROTEIN-NAMES", "reactome_term", "KINASE", 
                                    "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT", 
                                    "ON_OTHER_INTERACT", "REG_SITES_NOTES"))
    
    print("Step 1")
    step_1 <- scores_list$pwms_scores[[list_position]]$peptide_p |>
      pivot_longer(cols = !(contains("annotation") | contains("peptide")),
                   names_to ="kinase_gene",
                   values_to = "motif.p.value") |>
      left_join(motif_score_table, 
                by = c("annotation" = "annotation",
                       "peptide" = "peptide",
                       "kinase_gene" = "kinase_gene")) |>
      inner_join(up_or_down_kinases, 
                 by = c("kinase_gene" = "kinase_gene")) |>
      dplyr::filter(motif.p.value < 0.2 )
    
    print("Step 2")
    step_2 <- step_1 |>
      left_join(ks_tbl |>
                  mutate(KINASE = toupper(KINASE)) |>
                  dplyr::distinct( KINASE, KIN_ACC_ID),
                by = c("kinase_gene" = "KINASE")) |>
      left_join(ks_tbl |>
                  mutate(GENE = toupper(GENE)) |>
                  dplyr::distinct( GENE, KIN_ACC_ID),
                by = c("kinase_gene" = "GENE")) |>
      mutate(kinase_uniprot_acc = ifelse(is.na( KIN_ACC_ID.x),
                                         KIN_ACC_ID.y,
                                         KIN_ACC_ID.x)) |>
      dplyr::select(-KIN_ACC_ID.x, -KIN_ACC_ID.y) |>
      inner_join(fc_pval_tab |>
                   dplyr::filter(pval < p_value_cutoff),
                 by = c("annotation" = "annotation",
                        "peptide" = "peptide")) |>
      left_join(annotated_data |>
                  dplyr::filter(comparison == swing_out_list$comparison[[list_position]]) |>
                  dplyr::select(annotation, sites_id, peptide),
                by = c("annotation" = "annotation",
                       "peptide" = "peptide"))
    rm(step_1)
    gc()
    
    print("Step 3")
    step_3 <- step_2 |>
      left_join(phosphositeplus |>
                  distinct(kinase_gene, kinase_class), 
                by = c("kinase_gene" = "kinase_gene"))
    rm(step_2)
    gc()
    
    print("Step 4")
    plan(multisession, workers = num_cores)
    step_4 <- step_3 |>
      dplyr::mutate(substrate_gene_name = furrr::future_map(annotation,
                                                            \(x){ str_split(x, "\\|") |> 
                                                                purrr::map_chr(2)}))
    rm(step_3)
    gc()
    
    print("Step 5")
    step_5 <- step_4 |>
      left_join(de_phos |>
                  dplyr::filter(comparison == swing_out_list$comparison[[list_position]]) |>
                  dplyr::select(one_of(selected_columns)),
                by = c("sites_id" = "sites_id"))
    rm(step_4)
    gc()
    
    print("Step 6")
    selected_scores_list_help <- step_5 |>
      left_join(uniprot_kinases |>
                  dplyr::select(-one_of(c("KEYWORDS", "Keywords"))),
                by = c("kinase_gene" = "gene_name")) |>
      dplyr::rename(kinase_gene_name = "kinase_gene") |>
      distinct()
    rm(step_5)
    gc()
    
    if("KINASE" %in% colnames(selected_scores_list_help)) {
      selected_scores_list_help <- selected_scores_list_help |>
        dplyr::rename(one_known_kinase = "KINASE")
      print("Update KINASE")
    }
    
    if(!"one_known_kinase" %in% colnames(selected_scores_list_help)) {
      stop("Error, problem, column one_known_kinase is missing.")
    }
    
    print("Step 7")
    selected_scores_list_help <- selected_scores_list_help |>
      dplyr::mutate(known_upstream_kinase = one_known_kinase ) |>
      separate_rows(one_known_kinase , sep= "//") |>
      dplyr::mutate(known_upstream_kinase = toupper(known_upstream_kinase)) |>
      dplyr::mutate(prediction_match_known_kinase = case_when(kinase_class == known_upstream_kinase ~ TRUE,
                                                              TRUE ~ FALSE))
    
    print("Step 8")
    #- Use mouse or human uniprot accession if it makes sense to do so.
    selected_scores_list <- selected_scores_list_help
    if(length(intersect(is_phosphoylated_tbl$uniprot_acc, uniprot_kinases$uniprot_acc_human)) > 0) {
      print("human uniprot accession used")
      selected_scores_list <- selected_scores_list_help |>
        left_join(is_phosphoylated_tbl, 
                  by = c("uniprot_acc_human" = "uniprot_acc")) |>
        distinct()
    } else if(length(intersect(is_phosphoylated_tbl$uniprot_acc, uniprot_kinases$uniprot_acc_mouse) > 0 )) {
      print("mouse uniprot accession used")
      selected_scores_list <- selected_scores_list_help |>
        left_join(is_phosphoylated_tbl, 
                  by = c("uniprot_acc_mouse" = "uniprot_acc")) |>
        distinct()
    } else {
      print("neither human or mouse uniprot accession used")
      list_of_kinswinger_columns <- setdiff(list_of_kinswinger_columns,
                                            c("is_kinase_phosphorylated", "kinase_family", 
                                              "kinase_uniprot_id_human", "kinase_uniprot_acc_human", 
                                              "kinase_uniprot_id_mouse", "kinase_uniprot_acc_mouse", 
                                              "known_upstream_kinase", "prediction_match_known_kinase", 
                                              "reactome_term", "ON_FUNCTION", "ON_PROCESS",
                                              "ON_PROT_INTERACT", "ON_OTHER_INTERACT","REG_SITES_NOTES"))
    }
    
    print("Step 9")
    selected_scores_list_cln_step_1 <- selected_scores_list |>
      separate(annotation, 
               into = c("substrate_uniprot_acc", "subsrate_gene_symbol", 
                        "phosphosite_position", "sequence_context"), 
               sep = "\\|") |>
      dplyr::mutate(p_value_kinswingr = case_when(swing > 0 ~ p_greater,
                                                  swing < 0 ~ p_less,
                                                  TRUE ~ 1)) |>
      dplyr::rename(kinase_gene_symbol = "one_known_kinase",
                    substrate_name = "PROTEIN-NAMES",
                    phosphosite_log2FC = "fc",
                    phosphosite_fdr_value = "pval",
                    swing_kinswingr = "swing",
                    pos_kinswingr = "pos",
                    neg_kinswingr = "neg",
                    all_kinswingr = "all",
                    pk_kinswingr = "pk",
                    nk_kinswingr = "nk",
                    swing_raw_kinswingr = "swing_raw",
                    n_kinswingr = "n",
                    p_greater_kinswingr = "p_greater",
                    p_less_kinswingr = "p_less")
    
    selected_scores_list_cln_step_2 <- selected_scores_list_cln_step_1
    
    print("Step 10")
    if("Family" %in% colnames(selected_scores_list_cln_step_1)) {
      print("filter for human or mouse")
      selected_scores_list_cln_step_2 <- selected_scores_list_cln_step_1 |>
        dplyr::rename(kinase_family = "Family",
                      kinase_uniprot_id_human = "uniprot_id_human",
                      kinase_uniprot_acc_human = "uniprot_acc_human",
                      kinase_uniprot_id_mouse = "uniprot_id_mouse",
                      kinase_uniprot_acc_mouse = "uniprot_acc_mouse") |>
        dplyr::filter((taxonomy_id == 9606 &  kinase_uniprot_acc_human == kinase_uniprot_acc) | ## human only
                        (taxonomy_id == 10090 &  kinase_uniprot_acc_mouse == kinase_uniprot_acc) |  ## mouse only
                        !(taxonomy_id %in% c(9606, 10090)))
      
    }
    
    print("Step 11")
    if(file.exists(uniprot_to_gene_symbol_file)) {
      print("Add uniprot gene symbol")
      selected_scores_list_cln_step_3 <- selected_scores_list_cln_step_2 |>
        left_join(uniprot_tab_delimited_tbl |>
                    dplyr::select(!!rlang::sym(protein_id_lookup_column),
                                  !!rlang::sym(gene_symbol_column)) |>
                    dplyr::rename(kinsae_gene_name_uniprot = gene_symbol_column),
                  by = join_by(kinase_uniprot_acc == !!rlang::sym(protein_id_lookup_column)))
    } else {
      selected_scores_list_cln_step_3 <- selected_scores_list_cln_step_2
    }
    
    print("Step 12")
    print(colnames(selected_scores_list_cln_step_3))
    print(list_of_kinswinger_columns)
    selected_scores_list_cln_final <- selected_scores_list_cln_step_3 |>
      dplyr::select(any_of(list_of_kinswinger_columns))
    
    if(file.exists(uniprot_to_gene_symbol_file)) {
      selected_scores_list_cln_final <- selected_scores_list_cln_final |>
        relocate(kinsae_gene_name_uniprot, 
                 .before = "kinase_gene_symbol")
    }
    return(selected_scores_list_cln_final)
  }
  
  
  gc()
  
  selected_scores_list <- purrr::map(seq_along(swing_out_list$comparison),
                                     \(x) compileKinswingerResults(x))
  
  names(selected_scores_list) <- swing_out_list$comparison
  
  purrr::walk2(selected_scores_list,
               swing_out_list$comparison,
               \(.x, .y) { 
                 vroom::vroom_write(.x |>
                                      mutate(comparison = .y) |> 
                                      relocate(comparison, .before="substrate_uniprot_acc"), 
                                    file.path(output_dir,
                                              paste0("selected_kinase_substrate_",
                                                     kinase_specificity, 
                                                     "_", .y, ".tsv")))
               })
}