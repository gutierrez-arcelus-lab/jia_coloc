library(tidyverse)
library(readxl)

gwas <- read_excel("./TRAF1_region_summary_stats.xlsx")

gwas %>%
    select(start = pos_build36, snp = SNP) %>%
    mutate(chr = "chr9", end = start, start = start - 1) %>%
    select(chr, start, end, snp) %>%
    write_tsv("./data/gwas/TRAF1_region_hg18.bed", col_names = FALSE)
