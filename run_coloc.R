library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
import_eqtl <- safely(function(ftp_path, region, gene_ids) {
  
    #Fetch summary statistics with seqminer
    fetch_table <- 
	tabix.read.table(tabixFile = ftp_path, 
			 tabixRange = region) %>%
	as_tibble()

    column_names <- read_tsv(ftp_path, n_max = 1) %>%
	colnames()

    colnames(fetch_table) <- column_names

    fetch_table %>%
	filter(gene_id %in% gene_ids, type == "SNP") %>%
	filter(rsid != "NA\r") %>%
	add_count(gene_id, rsid) %>%
	filter(n == 1) %>%
	select(-n) %>%
	add_count(gene_id, position) %>%
	filter(n == 1) %>%
	select(-n) %>%
	mutate(rsid = sub("\\\r$", "", rsid))
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
paths2_df <- read_tsv("./data/coloc_inputs/eqtl_catalogue_paths.tsv")



# GWAS summ stats
gwas_data <- read_tsv("./data/coloc_inputs/gwas_data.tsv")

region <- read_lines("./data/coloc_inputs/region.txt")

# Genes in the region
genes_df <- read_tsv("./data/annotation/gene_annots.tsv")

# Import eQTL catalog
plan(multisession, workers = 8)

eqtl_database <- paths_df %>%
    mutate(data = future_map(ftp_path, import_eqtl, region = region, gene_ids = genes_df$gene_id))

error_df <- eqtl_database %>%
    mutate(error = map(data, "error")) %>%
    select(-ftp_path, -data)

write_rds(error_df, "./results/eqtl_catalog_query_log.rds")

success <- unlist(map(error_df$error, is.null))

eqtl_cat_df <- eqtl_database %>%
    filter(success) %>%
    mutate(data = map(data, "result")) %>%
    unnest(cols = c(data)) %>%
    mutate(eqtl_sample_size = an/2L,
	   eqtl_varbeta = se^2) %>%
    select(study, tissue = tissue_label, condition = condition_label, qtl_group, feature = quant_method, gene_id,
	   rsid, position, eqtl_sample_size, eqtl_maf = maf, eqtl_beta = beta, eqtl_varbeta)

write_tsv(eqtl_cat_df, "./results/eqtl_catalogue_retrieved_data.tsv.gz")

min_df <- inner_join(eqtl_cat_df, gwas_data, by = c("rsid", "position")) %>%
    filter(!is.na(eqtl_varbeta))

coloc_results <- min_df %>%
    unite("id", c(study, tissue, condition, qtl_group, feature, gene_id), sep = ",") %>%
    split(.$id) %>%
    future_map(run_coloc)

coloc_results_df <- map_df(coloc_results, "results", .id = "id") %>%
    as_tibble() %>%
    separate(id, c("study", "tissue", "condition", "qtl_group", "feature", "gene_id"), sep = ",") %>%
    left_join(genes_df) %>%
    select(study, tissue, condition, qtl_group, gene_id, gene_name, feature, everything())

coloc_results_summary <- map(coloc_results, "summary") %>%
    map_dfr(~as_tibble(., rownames = "stat"), .id = "id") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) %>%
    separate(id, c("study", "tissue", "condition", "qtl_group", "feature", "gene_id"), sep = ",") %>%
    left_join(genes_df) %>%
    select(study, tissue, condition, qtl_group, gene_id, gene_name, feature, everything())
   
write_tsv(coloc_results_df, "./results/coloc_results_df.tsv")
write_tsv(coloc_results_summary, "./results/coloc_results_summary.tsv")

