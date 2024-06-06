# need to creat an .rproj file in the main directory
base_dir <- here::here()

#----- Variable RUV3 Number of Unwanted Factors
proteome_de_var_ruv_k <- 3
#----- Variable RUV3 Number of Unwanted Factors
phosphoproteome_de_var_ruv_k <- 3

#----- Variable Average Technical Replicate ID, column name in design matrix
proteome_var_average_replicates_id <- NA
#----- Variable Replicate Group ID, column name in design matrix
proteome_var_replicate_group_id <- NA

#----- Variable Control Genes Q-value Threshold
phosphoproteome_de_var_control_genes_q_val_thresh <- 0.05


#----- Variable NCBI Taxonomy ID, M.musculus=10090, H.sapien=9606
species_taxon <- 9606
proteome_annotate_var_taxonomy_id <- species_taxon
proteome_collate_var_taxonomy_id <- species_taxon
phosphoproteome_annotate_var_taxonomy_id <- species_taxon
phosphoproteome_kinswingr_var_taxonomy_id <- species_taxon


phosphoproteome_annotate_dir_output <- "results/phosphoproteomics/annot_phos"


#----- Input Protein Sequence Fasta File Path with UniProt Header
proteome_clean_file_fasta <- "data/rawdata/HomoSapiens20240304CanIso.fasta"
#----- Input Protein Sequence Fasta File Path with UniProt Header
phosphoproteome_clean_file_fasta <- "data/rawdata/HomoSapiens20240304CanIso.fasta"

#----- Input Raw Counts
proteome_annotate_file_raw_counts <- "data/rawdata/proteinGroups.txt"
#----- Input Protein Abundance Data File
proteome_clean_file_raw_counts <- "data/rawdata/proteinGroups.txt"


#---------- Proteome
#----- Output Directory Path for Proteome Temporary Files
proteome_dir_tmp <- "results/proteomics/cache"
#----- Input Design Matrix File
proteome_file_design_matrix <- "proteomics_design.tab"
#----- Variable Sample ID, column name in design matrix
proteome_var_sample_id <- "Sample_ID"

#----- Variable Group ID, column name in design matrix
proteome_var_group_id <- "group"

#----- Variable Before Technical Replicate IDs were Averaged and Merged
proteome_var_before_avg_design_matrix <- NA
#----- Input Contrasts File
proteome_file_contrasts <- "proteomics_contrast.tab"
#----- Variable Formula String
proteome_var_formula_string <- "~ 0 + group "
#----- Variable Plots Format
proteome_var_plots_format <- list("pdf","png")

#---------- Proteome Clean
#----- Output Directory Path for Proteome Clean
proteome_clean_dir_output <- "results/proteomics/clean_proteins"
#----- Output Cleaned Counts Data File
proteome_clean_file_counts <- "counts_table_cleaned.tab"
#----- Output Cleaned MaxQuant Accessions File Name
proteome_clean_file_accession <- "cleaned_accession_to_protein_group.tab"
#----- Output Protein Sequence and Header Fasta File
proteome_clean_file_fasta_meta <- "aa_seq_tbl.RDS"

#----- Variable String Pattern to Identify Abundance Columns
proteome_clean_var_column_pattern <- "Reporter intensity corrected"
#----- Variable Regular Expression to Identify Group Abundance Columns
proteome_clean_var_group_pattern <- ""
#----- Variable Regular Expression to Identify Abundance Columns
proteome_clean_var_pattern_suffix <- "_\\d+"
#----- Variable Regular Expression to Identify Abundance Columns
proteome_clean_var_extract_suffix <- "_(\\d+)"
#----- Variable Razor and Unique Peptides Threshold per Group
proteome_clean_var_razor_unique_peptides_group_thresh <- 1
#----- Variable Unique Peptides Threshold per Group
proteome_clean_var_unique_peptides_group_thresh <- 1
#----- Variable Remove Reverse Decoy and Contaminant Peptides (logical)
proteome_clean_var_remove_more_peptides <- TRUE

#---------- Proteome DE
#----- Output Directory Path for Proteome DE
proteome_de_dir_output <- "results/proteomics/de_proteins"
#----- Variable Limma treat fuction for minimum log2 fold change
proteome_de_var_treat_lfc_cutoff <- NA
#----- Variable Normalization method, ["scale", "quantile" or "cyclicloess"]
proteome_de_var_normalization <- "cyclicloess"
#----- Variable Handling Zeros, ["mark_as_na", "add_pseudo_count"]
proteome_de_var_handling_zeros <- "mark_as_na"
#----- Variable Implement Imputation, [TRUE, FALSE]
proteome_de_var_imputation <- FALSE
#----- Variable Impute Minimum Percentage, 0 to 1 indicating percentage of quantified values
proteome_de_var_impute_min_percent <- NA
#----- Variable Impute Minimum Number of Groups, integer of replicates for percent filter 
proteome_de_var_impute_min_num_of_groups <- NA
#----- Variable Impute Specific Percentage, 0 to 1 indicating minimum percentage
proteome_de_var_impute_specific_percent <- NA
#----- Variable Remove Imputed, [TRUE, FALSE]
proteome_de_var_remove_imputed <- FALSE
#----- Variable Remove Percentage, 0 to 1 indicating minimum percentage
proteome_de_var_remove_imputed_perc <- 1
#----- Variable Empirical Bayes, [TRUE, FALSE]
proteome_de_var_eBayes_trend <- TRUE
#----- Variable Empirical Bayes Robust, [TRUE, FALSE]
proteome_de_var_eBayes_robust <- TRUE
#----- Variable Maximum Number of Samples with Missing Values per Group to Remove Protein
proteome_de_var_max_num_samples_miss_per_group <- 0
#----- Variable Abundance threshold
proteome_de_var_abundance_threshold <- 0
#----- Variable Regular Expression to Identify Group Abundance Columns
proteome_de_var_group_pattern <- "\\d+"
#----- Variable Q-value Threshold
proteome_de_var_q_val_thresh <- 0.05
#----- Variable Control Genes Q-value Threshold
proteome_de_var_control_genes_q_val_thresh <- 0.05

#----- Variable Number of Negative Control Proteins
proteome_de_var_num_neg_ctrl <- 100
#----- Variable RUV Method
proteome_de_var_ruv_method <- "ruv3"
#----- Input Protein Abundance Counts Table
proteome_de_file_counts_table <- "results/proteomics/clean_proteins/counts_table_cleaned.tab"
#----- Variable Row ID
proteome_de_var_row_id <- "uniprot_acc"
#----- Variable File Prefix
proteome_de_var_prefix <- "de_proteins"

#---------- Proteome Annotate
#----- Output Directory Path for Proteome Annotate
proteome_annotate_dir_output <- "results/proteomics/annot_proteins"

#----- Input Results Table Wide
proteome_annotate_file_input_wide <- "results/proteomics/de_proteins/de_proteins_wide.tsv"
#----- Input Results Table Long
proteome_annotate_file_input_long <- "results/proteomics/de_proteins/de_proteins_long.tsv"
#----- Input Protein Cleaned Accessions
proteome_annotate_file_ids <- "results/proteomics/clean_proteins/cleaned_accession_to_protein_group.tab"

#----- Output Results Table Wide
proteome_annotate_file_output_wide <- "de_proteins_wide_annot.tsv"
#----- Output Results Table Long
proteome_annotate_file_output_long <- "de_proteins_long_annot.tsv"
#----- Input Reactome File
proteome_annotate_file_reactome <- "UniProt2Reactome.txt"
#----- Input Uniprot File
proteome_annotate_file_uniprot <- "uniprot_data.RDS"

#---------- Proteome Graphs
#----- Output Directory Path for Proteome Annotate
proteome_graphs_dir_output <- "results/proteomics/publication_graphs"
#----- Input Directory Path for Proteome DE
proteome_graphs_dir_input <- "results/proteomics/de_proteins"
#----- Variable Row ID
proteome_graphs_var_row_id <- "uniprot_acc"
#----- Variable Q-value Threshold
proteome_graphs_var_q_val_thresh <- 0.05
#----- Input Proteome Annotate Results Table Long
proteome_graphs_file_de_proteins_long_file <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
#----- Variable Top Number of Genes
proteome_graphs_var_top_x_gene_name <- 5
#----- Variable Data Type, ["proteomics", "phosphoproteomics"]
proteome_graphs_var_data_type <- "proteomics"
#----- Variable log2FC Column Name
proteome_graphs_var_log2fc_column <- "log2FC"
#----- Variable FDR Column Name
proteome_graphs_var_fdr_column <- "q.mod"

#---------- Proteome Pathways
#----- Input Differentially expressed proteins
proteome_pathways_file_proteins <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
#----- Variable log2FC Column Name
proteome_pathways_var_log2fc_column <- "log2FC"
#----- Variable FDR Column Name
proteome_pathways_var_fdr_column <- "q.mod"
#----- Variable Q-value Threshold
proteome_pathways_var_q_val_thresh <- 0.05
#----- Variable Row ID
proteome_pathways_var_row_id <- "uniprot_acc"
#----- Variable P-value Threshold
proteome_pathways_var_p_val_thresh <- 0.05
#----- Input Uniprot Accession to Gene Symbol
proteome_pathways_file_uniprot_to_gene_symbol <- "data/uniprot/data.tab"
#----- Variable Uniprot Accession to Gene Symbol Protein ID Column
proteome_pathways_var_protein_id_column <- "Entry"
#----- Variable Uniprot Accession to Gene Symbol Protein Gene Names Column
proteome_pathways_var_gene_names_column <- "Gene Names"
#----- Variable Max Gene Set Size
proteome_pathways_var_max_gene_set_size <- list("250,500")
#----- Variable Min Gene Set Size
proteome_pathways_var_min_gene_set_size <- list("5,10")

#---------- Proteome Pathways GO
#----- Output Directory Path for Proteome Pathways
proteome_pathways_go_dir_output <- "results/proteomics/de_proteins_go_list"
#----- Input Annotation File
proteome_pathways_go_file_annotation <- "data/uniprot/go_terms_table_python_all.tab"
#----- Input Dictionary File
proteome_pathways_go_file_dictionary <- "data/uniprot/go_terms_table_python_all.tab"
#----- Variable Annotation Type, added to output file
proteome_pathways_go_var_annotation_type <- "go_enrichment"
#----- Variable Annotation ID, in annotation file
proteome_pathways_go_var_annotation_id <- "go_id"
#----- Variable Annotation Column, column name
proteome_pathways_go_var_annotation_column <- "go_term"
#----- Variable Aspect Column, column name
proteome_pathways_go_var_aspect_column <- "go_type"

#---------- Proteome Pathways Reactome
#----- Output Directory Path for Proteome Pathways
proteome_pathways_reactome_dir_output <- "results/proteomics/de_proteins_reactome_list"
#----- Input Annotation File
proteome_pathways_reactome_file_annotation <- "data/reactome/reactome_data_with_header.txt"
#----- Input Dictionary File
proteome_pathways_reactome_file_dictionary <- "data/reactome/reactome_data_with_header.txt"
#----- Variable Annotation Type, added to output file
proteome_pathways_reactome_var_annotation_type <- "reactome_enrichment"
#----- Variable Annotation ID, in annotation file
proteome_pathways_reactome_var_annotation_id <- "reactome_id"
#----- Variable Annotation Column, column name
proteome_pathways_reactome_var_annotation_column <- "reactome_term"
#----- Variable Aspect Column, column name
proteome_pathways_reactome_var_aspect_column <- NULL


#---------- Proteome Pathways Kegg
#----- Output Directory Path for Proteome Pathways
proteome_pathways_kegg_dir_output <- "results/proteomics/de_proteins_kegg_list"
#----- Input Annotation File
proteome_pathways_kegg_file_annotation <- "data/kegg/gene_sets.tab"
#----- Input Dictionary File
proteome_pathways_kegg_file_dictionary <- "data/kegg/gene_sets.tab"
#----- Variable Annotation Type, added to output file
proteome_pathways_kegg_var_annotation_type <- "kegg_enrichment"
#----- Variable Annotation ID, in annotation file
proteome_pathways_kegg_var_annotation_id <- "pathway_id"
#----- Variable Annotation Column, column name
proteome_pathways_kegg_var_annotation_column <- "pathway_name"
#----- Variable Aspect Column, column name
proteome_pathways_kegg_var_aspect_column <- NULL

#---------- Proteome Collate ----
## Proteome Collate GO ----
#----- Output Directory Path for Proteome Collate
proteome_collate_go_dir_output <- file.path(base_dir, "results/proteomics/Protein_Enrichment_Comparisons/GO")
  
#----- Input Enrichment File
proteome_collate_go_file_enrichment <- file.path(base_dir, "results/proteomics/de_proteins_go_list/*/go_enrichment_table_*.tab")
#----- Variable Minimum Number of Query Proteins
proteome_collate_go_var_min_num_query_proteins_matching <- 3
#----- Variable Selected Min Gene Set Size
proteome_collate_go_var_selected_min_set_size <- 5
#----- Variable Selected Max Gene Set Size
proteome_collate_go_var_selected_max_set_size <- 250
#----- Variable Run Revigo, [TRUE, FALSE]
proteome_collate_go_var_run_revigo <- TRUE
#----- Variable Input Type
proteome_collate_go_var_input_type <- NA
#----- Variable Revigo Cutoff
proteome_collate_go_var_revigo_cutoff <- 0.7
#----- Variable Annotation Type
proteome_collate_go_var_annotation_type <- "go_enrichment"
#----- Variable Plot Width
proteome_collate_go_var_plot_width <- 10
#----- Variable Plot Height
proteome_collate_go_var_plot_height <- 18

## Proteome Collate Reactome ----
#----- Output Directory Path for Proteome Collate
proteome_collate_reactome_dir_output <- file.path(base_dir, "results/proteomics/Protein_Enrichment_Comparisons/Reactome")
#----- Input Enrichment File
proteome_collate_reactome_file_enrichment <- file.path(base_dir, "results/proteomics/de_proteins_reactome_list/*/reactome_enrichment_table_*.tab")
#----- Variable Minimum Number of Query Proteins
proteome_collate_reactome_var_min_num_query_proteins_matching <- 3
#----- Variable Selected Min Gene Set Size
proteome_collate_reactome_var_selected_min_set_size <- 5
#----- Variable Selected Max Gene Set Size
proteome_collate_reactome_var_selected_max_set_size <- 250
#----- Variable Run Revigo, [TRUE, FALSE]
proteome_collate_reactome_var_run_revigo <- FALSE
#----- Variable Input Type
proteome_collate_reactome_var_input_type <- "Reactome"
#----- Variable Revigo Cutoff
proteome_collate_reactome_var_revigo_cutoff <- 0.5
#----- Variable Annotation Type
proteome_collate_reactome_var_annotation_type <- "reactome_enrichment"
#----- Variable Plot Width
proteome_collate_reactome_var_plot_width <- 15
#----- Variable Plot Height
proteome_collate_reactome_var_plot_height <- 30

## Proteome Collate Kegg ----
#----- Output Directory Path for Proteome Collate
proteome_collate_reactome_dir_output <- file.path(base_dir, "results/proteomics/Protein_Enrichment_Comparisons/KEGG")
#----- Input Enrichment File
proteome_collate_reactome_file_enrichment <- file.path(base_dir, "results/proteomics/de_proteins_kegg_list/*/kegg_enrichment_table_*.tab")
#----- Variable Minimum Number of Query Proteins
proteome_collate_reactome_var_min_num_query_proteins_matching <- 3
#----- Variable Selected Min Gene Set Size
proteome_collate_reactome_var_selected_min_set_size <- 5
#----- Variable Selected Max Gene Set Size
proteome_collate_reactome_var_selected_max_set_size <- 250
#----- Variable Run Revigo, [TRUE, FALSE]
proteome_collate_reactome_var_run_revigo <- FALSE
#----- Variable Input Type
proteome_collate_reactome_var_input_type <- "KEGG"
#----- Variable Revigo Cutoff
proteome_collate_reactome_var_revigo_cutoff <- 0.5
#----- Variable Annotation Type
proteome_collate_reactome_var_annotation_type <- "kegg_enrichment"
#----- Variable Plot Width
proteome_collate_reactome_var_plot_width <- 15
#----- Variable Plot Height
proteome_collate_reactome_var_plot_height <- 30

#---------- Phosphoproteome ----
#----- Output Directory Path for Phosphoproteome Temporary Files
phosphoproteome_dir_tmp <- "results/phosphoproteomics/cache"

#---------- Phosphoproteome Clean
#----- Output Directory Path for Phosphoproteome Clean
phosphoproteome_clean_dir_output <- "results/phosphoproteomics/clean_phos"
#----- Output Protein Sequence and Header Fasta File
phosphoproteome_clean_file_fasta_meta <- "aa_seq_tbl.RDS"

#----- Input Phosphoprotein Abundance Data File
phosphoproteome_clean_file_raw_counts <- "data/rawdata/evidence.txt"
#----- Variable String Pattern to Identify Abundance Columns
phosphoproteome_clean_var_column_pattern <- "Reporter intensity corrected"
#----- Variable List of Additional Columns
phosphoproteome_clean_var_add_cols_string <- ""
#----- Variable Regular Expression to Identify Abundance Columns
phosphoproteome_clean_var_pattern_suffix <- "_\\d+"
#----- Variable Regular Expression to Identify Abundance Columns
phosphoproteome_clean_var_extract_suffix <- "_(\\d+)"
#----- Variable Site Probability Threshold
phosphoproteome_clean_var_site_prob_threshold <- 0.75
#----- Variable Recover Site Probability Threshold
phosphoproteome_clean_var_recover_site_prob_thresh <- 0.5

#---------- Phosphoproteome DE
#----- Output Directory Path for Phosphoproteome DE
phosphoproteome_de_dir_output <- "results/phosphoproteomics/de_phos"
#----- Variable Limma treat fuction for minimum log2 fold change
phosphoproteome_de_var_treat_lfc_cutoff <- NA
#----- Variable Normalization method, ["scale", "quantile" or "cyclicloess"]
phosphoproteome_de_var_normalization <- "cyclicloess"
#----- Variable Handling Zeros, ["mark_as_na", "add_pseudo_count"]
phosphoproteome_de_var_handling_zeros <- "mark_as_na"
#----- Variable Implement Imputation, [TRUE, FALSE]
phosphoproteome_de_var_imputation <- FALSE
#----- Variable Impute Minimum Percentage, 0 to 1 indicating percentage of quantified values
phosphoproteome_de_var_impute_min_percent <- NA
#----- Variable Impute Minimum Number of Groups, integer of replicates for percent filter 
phosphoproteome_de_var_impute_min_num_of_groups <- NA
#----- Variable Impute Specific Percentage, 0 to 1 indicating minimum percentage
phosphoproteome_de_var_impute_specific_percent <- NA
#----- Variable Remove Imputed, [TRUE, FALSE]
phosphoproteome_de_var_remove_imputed <- FALSE
#----- Variable Remove Percentage, 0 to 1 indicating minimum percentage
phosphoproteome_de_var_remove_imputed_perc <- 1
#----- Variable Empirical Bayes, [TRUE, FALSE]
phosphoproteome_de_var_eBayes_trend <- TRUE
#----- Variable Empirical Bayes Robust, [TRUE, FALSE]
phosphoproteome_de_var_eBayes_robust <- TRUE
#----- Variable Maximum Number of Samples with Missing Values per Group to Remove Protein
phosphoproteome_de_var_max_num_samples_miss_per_group <- 0
#----- Variable Abundance threshold
phosphoproteome_de_var_abundance_threshold <- 0
#----- Variable Regular Expression to Identify Group Abundance Columns
phosphoproteome_de_var_group_pattern <- "\\d+"
#----- Variable Q-value Threshold
phosphoproteome_de_var_q_val_thresh <- 0.05

#----- Variable Number of Negative Control Proteins
phosphoproteome_de_var_num_neg_ctrl <- 100
#----- Variable RUV Method
phosphoproteome_de_var_ruv_method <- "ruv3"
#----- Input Protein Abundance Counts Table
phosphoproteome_de_file_counts_table <- "results/phosphoproteomics/clean_phos/sum_phoshpsites.tsv"
#----- Variable Row ID
phosphoproteome_de_var_row_id <- "sites_id"
#----- Variable File Prefix
phosphoproteome_de_var_prefix <- "de_phos"

#---------- Phosphoproteome Annotate
#----- Output Directory Path for Phosphoproteome Annotate

#----- Input Results Table Wide
phosphoproteome_annotate_file_input_wide <- "results/phosphoproteomics/de_phos/de_phos_wide.tsv"
#----- Input Results Table Long
phosphoproteome_annotate_file_input_long <- "results/phosphoproteomics/de_phos/de_phos_long.tsv"
#----- Input PhosphoSitePlus Database
phosphoproteome_annotate_dir_phosphosite_db <- "data/phosphosites"
#----- Input Raw Counts
phosphoproteome_annotate_file_raw_counts <- "results/phosphoproteomics/clean_phos/sum_phoshpsites.tsv"
#----- Output Results Table Wide
phosphoproteome_annotate_file_output_wide <- "de_phos_wide_annot.tsv"
#----- Output Results Table Long
phosphoproteome_annotate_file_output_long <- "de_phos_long_annot.tsv"
#----- Variable Near PTM Number of Residues
phosphoproteome_annotate_var_near_ptm_num_residues <- 5
#----- Input Reactome File
phosphoproteome_annotate_file_reactome <- "reactome_data.txt"
#----- Input Uniprot File
phosphoproteome_annotate_file_uniprot <- "uniprot_data_phos.RDS"

#---------- Phosphoproteome Normalize
#----- Output Directory Path for Phosphoproteome Normalize
phosphoproteome_normalize_dir_output <- "results/phosphoproteomics/norm_phos_by_prot_abundance"
#----- Input Differentially Expressed Proteins
phosphoproteome_normalize_file_proteins <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
#----- Input Differentially Expressed Phosphosites
phosphoproteome_normalize_file_phospho <- "results/phosphoproteomics/annot_phos/de_phos_long_annot.tsv"
#----- Variable Protein Group Filter, ["protein_rank", "best_p_value" or "best_log_fc"]
phosphoproteome_normalize_var_protein_group_filter <- "protein_rank"

#---------- Phosphoproteome Graphs
#----- Output Directory Path for Proteome Annotate
phosphoproteome_graphs_dir_output <- "results/phosphoproteomics/publication_graphs_phos"
#----- Input Directory Path for Proteome DE
phosphoproteome_graphs_dir_input <- "results/phosphoproteomics/de_phos"
#----- Variable Row ID
phosphoproteome_graphs_var_row_id <- "sites_id"
#----- Variable Q-value Threshold
phosphoproteome_graphs_var_q_val_thresh <- 0.05
#----- Input Proteome Annotate Results Table Long
phosphoproteome_graphs_file_de_proteins_long_file <- "results/phosphoproteomics/annot_phos/de_phos_long_annot.tsv"
#----- Variable Top Number of Genes
phosphoproteome_graphs_var_top_x_gene_name <- 5
#----- Variable Data Type, ["proteomics", "phosphoproteomics"]
phosphoproteome_graphs_var_data_type <- "phosphoproteomics"
#----- Variable log2FC Column Name
phosphoproteome_graphs_var_log2fc_column <- "log2FC"
#----- Variable FDR Column Name
phosphoproteome_graphs_var_fdr_column <- "q.mod"

#---------- Phosphoproteome Pathways
#----- Input Differentially expressed proteins
phosphoproteome_pathways_file_proteins <- "results/proteomics/annot_proteins/de_proteins_long_annot.tsv"
#----- Variable log2FC Column Name
phosphoproteome_pathways_var_log2fc_column <- "norm_phos_logFC"
#----- Variable FDR Column Name
phosphoproteome_pathways_var_fdr_column <- "combined_q_mod"
#----- Variable Row ID
phosphoproteome_pathways_var_row_id <- "uniprot_acc"
#----- Variable P-value Threshold
phosphoproteome_pathways_var_p_val_thresh <- 0.05
#----- Variable Site P-value Threshold
phosphoproteome_pathways_var_site_p_val_thresh <- 0.05
#----- Input Uniprot Accession to Gene Symbol
phosphoproteome_pathways_file_uniprot_to_gene_symbol <- "data/uniprot/data.tab"
#----- Variable Uniprot Accession to Gene Symbol Protein ID Column
phosphoproteome_pathways_var_protein_id_column <- "Entry"
#----- Variable Uniprot Accession to Gene Symbol Protein Gene Names Column
phosphoproteome_pathways_var_gene_names_column <- "Gene Names"
#----- Variable Max Gene Set Size
phosphoproteome_pathways_var_max_gene_set_size <- list("250,500")
#----- Variable Min Gene Set Size
phosphoproteome_pathways_var_min_gene_set_size <- list("5,10")
#----- Variable Background, ["phosphoproteins" or "proteins_and_phosphoproteins"]
phosphoproteome_pathways_var_background <- "proteins_and_phosphoproteins"
#----- Input Phosphorylation logFC Normalized by Protein logFC
phosphoproteome_pathways_file_norm_phos_logfc_file <- "results/phosphoproteomics/norm_phos_by_prot_abundance/norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv"
#----- Variable Gene Name Column in Results File
phosphoproteome_pathways_var_results_gene_name_column <- "gene_name"

#---------- Phosphoproteome Pathways GO
#----- Output Directory Path for Proteome Pathways
phosphoproteome_pathways_go_dir_output <-  "results/phosphoproteomics/phosphoproteins_go"
#----- Input Annotation File
phosphoproteome_pathways_go_file_annotation <- "data/uniprot/go_terms_table_python_all.tab"
#----- Input Dictionary File
phosphoproteome_pathways_go_file_dictionary <- "data/uniprot/go_terms_table_python_all.tab"
#----- Variable Annotation Type, added to output file
phosphoproteome_pathways_go_var_annotation_type <- "go_enrichment"
#----- Variable Annotation ID, in annotation file
phosphoproteome_pathways_go_var_annotation_id <- "go_id"
#----- Variable Annotation Column, column name
phosphoproteome_pathways_go_var_annotation_column <- "go_term"
#----- Variable Aspect Column, column name
phosphoproteome_pathways_go_var_aspect_column <- "go_type"

#---------- Phosphoproteome Pathways Reactome
#----- Output Directory Path for Proteome Pathways
phosphoproteome_pathways_reactome_dir_output <- "results/phosphoproteomics/phosphoproteins_reactome"
#----- Input Annotation File
phosphoproteome_pathways_reactome_file_annotation <- "data/reactome/reactome_data_with_header.txt"
#----- Input Dictionary File
phosphoproteome_pathways_reactome_file_dictionary <- "data/reactome/reactome_data_with_header.txt"
#----- Variable Annotation Type, added to output file
phosphoproteome_pathways_reactome_var_annotation_type <- "reactome_enrichment"
#----- Variable Annotation ID, in annotation file
phosphoproteome_pathways_reactome_var_annotation_id <- "reactome_id"
#----- Variable Annotation Column, column name
phosphoproteome_pathways_reactome_var_annotation_column <- "reactome_term"
#----- Variable Aspect Column, column name
phosphoproteome_pathways_reactome_var_aspect_column <- NULL

#---------- Phosphoproteome Pathways Kegg
#----- Output Directory Path for Proteome Pathways
phosphoproteome_pathways_kegg_dir_output <- "results/phosphoproteomics/phosphoproteins_kegg"
#----- Input Annotation File
phosphoproteome_pathways_kegg_file_annotation <- "data/kegg/gene_sets.tab"
#----- Input Dictionary File
phosphoproteome_pathways_kegg_file_dictionary <- "data/kegg/gene_sets.tab"
#----- Variable Annotation Type, added to output file
phosphoproteome_pathways_kegg_var_annotation_type <- "kegg_enrichment"
#----- Variable Annotation ID, in annotation file
phosphoproteome_pathways_kegg_var_annotation_id <- "pathway_id"
#----- Variable Annotation Column, column name
phosphoproteome_pathways_kegg_var_annotation_column <- "pathway_name"
#----- Variable Aspect Column, column name
phosphoproteome_pathways_kegg_var_aspect_column <- NULL

#---------- Phosphoproteome Kinswingr
#----- Variable Number of Cores
phosphoproteome_kinswingr_var_num_cores <- 1
#----- Variable Random Seed
phosphoproteome_kinswingr_var_random_seed <- 123456
#----- Variable Number of Iteration for Scoring each Substrate against Kinase Motif
phosphoproteome_kinswingr_var_motif_score_iteration <- 1000
#----- Variable Number of Iteration for Swing Scoring
phosphoproteome_kinswingr_var_swing_iteration <- 1000
#----- Variable Minimum Number of Known Substrates for each Kinase
phosphoproteome_kinswingr_var_min_num_sites_per_kinase <- 10
#----- Input Uniprot Accession to Gene Symbol
phosphoproteome_kinswingr_file_uniprot_to_gene_symbol <- "data/uniprot/data.tab"
#----- Variable Uniprot Accession to Gene Symbol Protein ID Column
phosphoproteome_kinswingr_var_protein_id_column <- "Entry"
#----- Variable Uniprot Accession to Gene Symbol Protein Gene Names Column
phosphoproteome_kinswingr_var_gene_names_column <- "Gene Names"
#----- Variable P-value Cutoff for Kinase Significance
phosphoproteome_kinswingr_var_p_val_cutoff <- 0.05
#----- Variable P-value Cutoff for Kinase Volcano Plot
phosphoproteome_kinswingr_var_p_val_volcano <- 0.1
#----- Input PhosphoSitePlus Database
phosphoproteome_kinswingr_dir_phosphosite_db <- "data/phosphosites"
#----- Input Uniprot Kinase File
phosphoproteome_kinswingr_file_uniprot_kinase_file <- "data/uniprot/pkinfam.tab"
#----- Input Phosphorylation logFC Normalized by Protein logFC
phosphoproteome_kinswingr_file_norm_phos_logfc_file <- "results/phosphoproteomics/norm_phos_by_prot_abundance/norm_phosphosite_lfc_minus_protein_lfc_annotated.tsv"
#----- Output UniProt Accession of Atypical and Other Kinases
phosphoproteome_kinswingr_file_uniprot_other_kinase_file <- "atypical_and_other_kinases_data.RDS"
#----- Variable log2FC Column Name
phosphoproteome_kinswingr_var_log2fc_column <- "norm_phos_logFC"

#----- Variable FDR Column Name
phosphoproteome_kinswingr_var_fdr_column <- "combined_q_mod"
#----- Variable Include both single-site and multi-site phosphorylation 
#- Include both single-site and multi-site phosphorylation [NA], 
#- Include only single-site phosphorylation [FALSE]
#- Include only multi-site phosphorylation [TRUE]
phosphoproteome_kinswingr_var_single_multi_site <- NA
#----- Variable Reuse Previous Saved Result [TRUE or FALSE]
phosphoproteome_kinswingr_var_reuse_old <- FALSE
#----- Variable Plots Format
proteome_var_plots_format <- list("pdf","png")

#----- Phosphoproteome Kinswingr ST
#----- Output Directory Path for Phosphoproteome Kinswingr
phosphoproteome_kinswingr_st_dir_output <- "results/phosphoproteomics/phos_kinswingr_ST"
#----- Variable Kinase Specificity, Ser/Thr kinase = ST, Tyr kinase = Y, Ser/Thr/Tyr kinase = STY   
phosphoproteome_kinswingr_st_var_kinase_specificity <- "ST"

#----- Phosphoproteome Kinswingr STY
#----- Output Directory Path for Phosphoproteome Kinswingr
phosphoproteome_kinswingr_sty_dir_output <- "results/phosphoproteomics/phos_kinswingr_STY"
#----- Variable Kinase Specificity, Ser/Thr kinase = ST, Tyr kinase = Y, Ser/Thr/Tyr kinase = STY   
phosphoproteome_kinswingr_sty_var_kinase_specificity <- "STY"

#----- Phosphoproteome Kinswingr Y
#----- Output Directory Path for Phosphoproteome Kinswingr
phosphoproteome_kinswingr_y_dir_output <- "results/phosphoproteomics/phos_kinswingr_Y"
#----- Variable Kinase Specificity, Ser/Thr kinase = ST, Tyr kinase = Y, Ser/Thr/Tyr kinase = STY   
phosphoproteome_kinswingr_y_var_kinase_specificity <- "Y"




#---------- Variables Table
# knitr::kable(data.frame(Variable = ls(),
#                         Input = lapply(ls(), get) %>% unlist()),
#              caption = "Variables Table") %>%
#   kable_styling(bootstrap_options = c("striped","bordered", "hover",
#                                       "condensed", "responsive"))

# getVariables <- function(list_multi_variables) {
#   #- Vector of All variables
#   all_variables <- ls(envir=.GlobalEnv)[-match("getVariables", ls(envir=.GlobalEnv))]
#   #- Duplicate for Multi Variables
#   all_variables_extra <- all_variables
#   for(i in list_multi_variables) {
#     # print(paste0("Multivariable: ", i))
#     # print(paste0("Multivariable length: ", length(get(sym(i)))))
#     for(j in 2:length(get(sym(i)))) {
#       all_variables_extra <- all_variables_extra %>% 
#         append(i, match(i, all_variables_extra)[1])
#     }
#   }
#   all_variables_df <- data.frame(Variable = all_variables_extra, 
#                                  Input = lapply(all_variables, get) %>% unlist())
#   return(all_variables_df)
# }
# getVariables(list_multi_variables = c("proteome_var_plots_format",
#                                       "proteome_pathways_var_min_gene_set_size",
#                                       "proteome_pathways_var_max_gene_set_size",
#                                       "phosphoproteome_pathways_var_min_gene_set_size",
#                                       "phosphoproteome_pathways_var_max_gene_set_size")) %>% 
#   knitr::kable(caption = "Variables Table") %>%
#   kable_styling(bootstrap_options = c("striped","bordered", "hover",
#                                       "condensed", "responsive"))