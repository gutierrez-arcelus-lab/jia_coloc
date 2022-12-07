library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
import_eqtl <- safely(function(ftp_path, region, genes, header) {
  
    #Fetch summary statistics with seqminer
    fetch_table <- 
	tabix.read.table(tabixFile = ftp_path, 
			 tabixRange = region) %>%
	as_tibble()

    colnames(fetch_table) <- header
    
    filter(fetch_table, gene_id %in% genes)
})

run_coloc <- function(min_df) {

    eqtl_dataset <- list(beta = min_df$eqtl_beta,
			 varbeta = min_df$eqtl_varbeta,
			 N = min_df$eqtl_sample_size,
			 MAF = min_df$eqtl_maf,
			 type = "quant",
			 snp = min_df$rsid)

    gwas_dataset <- list(beta = min_df$gwas_beta,
			 varbeta = min_df$gwas_varbeta,
			 type = "cc",
			 snp = min_df$rsid)

    coloc_res <- coloc.abf(dataset1 = eqtl_dataset, 
			   dataset2 = gwas_dataset)

    coloc_res
}


# eQTL catalogue paths
paths_df <- read_tsv("./data/coloc_inputs/eqtl_catalogue_paths.tsv")

column_names <- read_tsv(paths_df$ftp_path[1], n_max = 1) %>%
    colnames()

region <- read_lines("./data/coloc_inputs/region.txt")

# Genes in the region
genes_df <- read_tsv("./data/coloc_inputs/gene_annots.tsv")

# Import eQTL catalog
plan(multisession, workers = 8)

eqtl_database <- paths_df %>%
    mutate(data = future_map(ftp_path, import_eqtl, 
			     region = region, genes = genes_df$gene_id, header = column_names))

# Failed data acquisition:
error_df <- eqtl_database %>%
    mutate(error = map(data, "error"))

success <- unlist(map(error_df$error, is.null))

write_rds(eqtl_database, "./results/eqtl_catalogue_retrieved_data.rds")

#eqtl_database <- read_rds("./results/eqtl_catalogue_retrieved_data.rds")

# Data filtering
eqtl_data_filter <- eqtl_database %>%
    select(-ftp_path) %>%
    mutate(data = map(data, "result")) %>%
    unnest(cols = c(data)) %>%
    filter(type == "SNP") %>%
    mutate(rsid = sub("\\\r$", "", rsid)) %>%
    filter(rsid != "NA") %>% 
    group_by(study, tissue_label, condition_label, qtl_group, molecular_trait_id, gene_id) %>%
    filter(any(pvalue < 5e-5)) %>%
    ungroup() %>%
    add_count(study, tissue_label, condition_label, qtl_group, molecular_trait_id, gene_id, rsid) %>%
    filter(n == 1) %>%
    select(-n) %>%
    add_count(study, tissue_label, condition_label, qtl_group, molecular_trait_id, gene_id, position) %>%
    filter(n == 1) %>%
    select(-n)

eqtl_cat_df <- eqtl_data_filter %>%
    mutate(eqtl_sample_size = an/2L,
	   eqtl_varbeta = se^2) %>%
    select(study, tissue = tissue_label, condition = condition_label, qtl_group, 
	   feature = quant_method, gene_id, molecular_trait_id,
	   rsid, position, eqtl_sample_size, eqtl_maf = maf, eqtl_beta = beta, eqtl_varbeta)

# GWAS summ stats
gwas_data <- read_tsv("./data/coloc_inputs/gwas_data.tsv")


# Run coloc
min_df <- inner_join(eqtl_cat_df, gwas_data, by = c("rsid", "position")) %>%
    filter(!is.na(eqtl_varbeta)) %>%
    filter(gwas_beta != 0, gwas_varbeta != 0)

coloc_results <- min_df %>%
    unite("id", c(study, tissue, condition, qtl_group, feature, gene_id, molecular_trait_id), sep = ",") %>%
    split(.$id) %>% 
    map(run_coloc)

coloc_results_df <- map_df(coloc_results, "results", .id = "id") %>%
    as_tibble() %>%
    separate(id, c("study", "tissue", "condition", "qtl_group", "feature", "gene_id", "molecular_trait_id"), sep = ",") %>%
    left_join(genes_df) %>%
    select(study, tissue, condition, qtl_group, gene_id, gene_name, feature, everything())

coloc_results_summary <- map(coloc_results, "summary") %>%
    map_dfr(~as_tibble(., rownames = "stat"), .id = "id") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) %>%
    separate(id, c("study", "tissue", "condition", "qtl_group", "feature", "gene_id", "molecular_trait_id"), sep = ",") %>%
    left_join(genes_df) %>%
    select(study, tissue, condition, qtl_group, gene_id, gene_name, feature, everything())

# Save results
write_tsv(coloc_results_df, "./results/coloc_results_df.tsv")
write_tsv(coloc_results_summary, "./results/coloc_results_summary.tsv")

