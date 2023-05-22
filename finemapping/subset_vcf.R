library(tidyverse)
library(readxl)

cmdargs <- commandArgs(TRUE)
vcfin <- cmdargs[1]
vcfout <- cmdargs[2]

work_dir <- "/lab-share/IM-Gutierrez-e2/Public/vitor/jia"

# GWAS summ stats
gwas_stats <- 
    file.path(work_dir, "data/gwas/TRAF1_region_summary_stats.xlsx") |>
    read_excel()

# read file after running liftover
gwas_stats_38 <- 
    file.path(work_dir, "data/gwas/TRAF1_region_hg38.bed") |>
    read_tsv(col_names = c("chr", "start", "end", "rsid")) |>
    select(chr, pos = end, rsid) |>
    left_join(gwas_stats, by = c("rsid" = "SNP"))

valid <- file.path(work_dir, "finemapping/data/valid_positions.txt") |>
    read_tsv()

# ref panel
Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)
tmp <- read_lines(vcfin, n_max = 5000)
header <- keep(tmp, grepl("^##", tmp))
panel <- read_tsv(vcfin, comment = "##")

# filter intersect
panel_filter <- panel |>
    inner_join(select(gwas_stats_38, chr, pos), 
	       join_by(`#CHROM` == chr, POS == pos)) |>
    inner_join(valid, join_by(`#CHROM` == chr, POS == pos))

# Write output
vcfout_tmp <- sub("\\.gz$", "", vcfout)
write_lines(header, vcfout_tmp)
write_tsv(panel_filter, vcfout_tmp, col_names = TRUE, append = TRUE)
unlink(vcfout)
unlink(paste0(vcfout, ".tbi"))
system(paste("bgzip", vcfout_tmp))
system(paste("tabix -p vcf", vcfout))
