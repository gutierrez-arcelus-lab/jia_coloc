library(tidyverse)
library(readxl)

# functions
download_tbi <- function(ftp) {
    
    ftptbi <- paste0(ftp, ".tbi")
    idx_file <- file.path(getwd(), basename(ftptbi))

    if (! file.exists(idx_file) ) download.file(ftptbi, dest = idx_file)
}


# eQTL catalogue tabix data
tabix_paths <- 
    "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv" %>%
    read_tsv()

imported_tabix_paths <- 
    "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv" %>%
    read_tsv()

paths_df <- bind_rows(select(tabix_paths, -sample_size), imported_tabix_paths)

paths_df %>%
    select(study, tissue_label, condition_label, qtl_group, quant_method, ftp_path) %>%
    write_tsv("./data/coloc_inputs/eqtl_catalogue_paths.tsv")

# download index files
walk(paths_df$ftp_path, download_tbi)


# GWAS summ stats
gwas_stats <- read_excel("./data/gwas/TRAF1_region_summary_stats.xlsx")

gwas_stats_38 <- "./data/gwas/TRAF1_region_hg38.bed" %>%
    read_tsv(col_names = c("chr", "start", "end", "rsid")) %>%
    select(chr, pos = end, rsid) %>%
    left_join(gwas_stats, by = c("rsid" = "SNP"))

topsnp <- gwas_stats_38 %>%
    arrange(additive_pvalue) %>%
    slice(1) %>%
    pull(pos)

coords <- gwas_stats_38 %>%
    filter(between(pos, topsnp - 5e5, topsnp + 5e5)) %>%
    summarise(region = range(pos)) %>%
    pull(region)

region <- sprintf("%s:%d-%d", sub("chr", "", gwas_stats_38$chr[1]), coords[1], coords[2])

write_lines(region, "./data/coloc_inputs/region.txt")

gwas_df <- gwas_stats_38 %>%
    filter(between(pos, coords[1], coords[2])) %>%
    select(rsid, pos, n_case, n_ctrl, or = additive_odds_ratio, pvalue = additive_pvalue) %>%
    mutate(log_or = log(or),
	   varbeta = (log_or^2) / qchisq(pvalue, df = 1, lower = FALSE),
	   type = "cc",
	   s = n_case/n_ctrl) %>%
    select(rsid, pos, log_or, varbeta, type, s)

gwas_min_df <- gwas_df %>%
    filter(log_or > 0) %>%
    select(rsid, position = pos, gwas_beta = log_or, gwas_varbeta = varbeta)

write_tsv(gwas_min_df, "./data/coloc_inputs/gwas_data.tsv")
