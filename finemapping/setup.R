library(tidyverse)
library(readxl)

# GWAS summ stats
gwas_stats <- read_excel("../data/gwas/TRAF1_region_summary_stats.xlsx")

# read file after running liftover
gwas_stats_38 <- "../data/gwas/TRAF1_region_hg38.bed" |>
    read_tsv(col_names = c("chr", "start", "end", "rsid")) |>
    select(chr, pos = end, rsid) |>
    left_join(gwas_stats, by = c("rsid" = "SNP"))

region <- range(gwas_stats_38$pos) |>
    paste(collapse = "-")

region <- paste(unique(gwas_stats_38$chr), region, sep = ":")

write_lines(region, "./data/gwas/gwas_region.txt")

# not sure if I can run susie with p-values from the dominance model, 
# even though sometimes it is the best model.
# see: https://github.com/stephenslab/susieR/issues/187
stats_df <- gwas_stats_38 |>
    select(chr, pos, rsid, major = major_allele, minor = minor_allele, 
	   n_case, n_ctrl, or = additive_odds_ratio, p = additive_pvalue) |>
    filter(or != 1) |>
    mutate(beta = log(or),
	   varbeta = (beta^2) / qchisq(p, df = 1, lower = FALSE),
	   se = sqrt(varbeta),
	   z = beta/se)

stats_df |>
    select(chr, pos) |>
    write_tsv("./data/valid_positions.txt")

