library(coloc)
library(GenomicRanges)
library(seqminer)
library(tidyverse)
library(furrr)

# functions
fetch <- function(f, r) {
    tabix.read.table(tabixFile = f, 
		     tabixRange = r) |>
    as_tibble()
}

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

main <- function(ftp, region, gwas_sum, genes) {

    header <- read_tsv(ftp, n_max = 1) |>
	names()
    
    eqtl <- fetch(f = ftp, r = region) |>
	as_tibble() |>
	setNames(header) |>
	filter(gene_id %in% genes)

    eqtl_filtered <- eqtl |>
	mutate(rsid = sub("\\\r$", "", rsid)) |>
	filter(rsid != "NA") |>
	group_by(molecular_trait_id, gene_id) |>
	filter(any(pvalue < 5e-5)) |>
	ungroup() |>
	add_count(molecular_trait_id, gene_id, rsid) |>
	filter(n == 1) |>
	select(-n) |>
	add_count(molecular_trait_id, gene_id, type, position) |>
	filter(n == 1) |>
	select(-n)

    if ( nrow(eqtl_filtered) == 0 ) {
	return(tibble(gene_id = NA, molecular_trait_id = NA, nsnps = NA, 
		      h0 = NA, h1 = NA, h2 = NA, h3 = NA, h4 = NA))
    }

    eqtl_cat_df <- eqtl_filtered |>
	mutate(eqtl_sample_size = an/2L,
	       eqtl_varbeta = se^2) |>
	select(gene_id, molecular_trait_id, rsid, position, 
	       eqtl_sample_size, eqtl_maf = maf, eqtl_beta = beta, eqtl_varbeta)

    min_df <- inner_join(eqtl_cat_df, gwas_sum, by = c("rsid", "position")) |>
	filter(!is.na(eqtl_varbeta)) |>
	filter(gwas_beta != 0, gwas_varbeta != 0)

    coloc_results <- min_df |>
	unite("id", c(gene_id, molecular_trait_id), sep = ",") |>
	{function(x) split(x, x$id)}() |>
	map(run_coloc)

    coloc_results_summary <- map(coloc_results, "summary") |>
	map_dfr(function(x) as_tibble(x, rownames = "stat"), .id = "id") |>
	pivot_wider(names_from = stat, values_from = value) |>
	select(id, nsnps, h0 = PP.H0.abf, h1 = PP.H1.abf, h2 = PP.H2.abf, h3 = PP.H3.abf, h4 = PP.H4.abf) |>
	separate(id, c("gene_id", "molecular_trait_id"), sep = ",") |>
	select(gene_id, molecular_trait_id, everything())

    coloc_results_summary
}


# eQTL catalogue URLs
ftps_df <- "data/coloc_inputs/eqtl_catalogue_paths.tsv" |>
    read_tsv() |>
    filter(quant_method == "ge", study != "GTEx")
    
# Region
region <- read_lines("./data/coloc_inputs/region.txt")

# genes to consider (TSS 250kb from GWAS hit)
genes_dat <- "data/coloc_inputs/gene_annots.tsv" |>
    read_tsv()

# GWAS summ stats
gwas_data <- "data/coloc_inputs/gwas_data.tsv" |>
    read_tsv()

# Set up parallel architecture
plan(cluster, workers = length(availableWorkers()))

# Run analysis
coloc_res_list <- list()
error_i <- seq_len(nrow(ftps_df))
i <- 1
    
while ( length(error_i) > 0 && i <= 30) {

    coloc_res_list[[i]] <- ftps_df |>
	slice(error_i) |> 
	mutate(res = future_map(ftp_path,  
				safely(function(x) main(ftp = x, 
							region = region, 
							gwas_sum = gwas_data,
							genes = genes_dat$gene_id))))

    errors <- coloc_res_list[[i]] |>
	mutate(err = map(res, "error")) |>
	pull(err)

    error_i <- errors |>
	map(~!is.null(.)) |>
	unlist() |>
	which()

    i <- i + 1
}

if ( length(error_i) > 0) stop("Could not retrieve or process eQTL catalogue data.")

results <- coloc_res_list |>
    bind_rows() |>
    mutate(dat = map(res, "result")) |>
    select(study, tissue_label, condition_label, qtl_group, quant_method, dat) |>
    unnest(cols = dat) |>
    left_join(genes_dat, by = "gene_id") |>
    select(study:quant_method, gene_id, gene_name, everything())

out_name <- "results/coloc_results_df_v2.tsv"

write_tsv(results, out_name)
