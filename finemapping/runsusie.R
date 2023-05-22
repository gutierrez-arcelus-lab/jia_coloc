library(tidyverse)
library(readxl)
library(susieR)

# GWAS summ stats
gwas_stats <- read_excel("./data/gwas/TRAF1_region_summary_stats.xlsx")

# read file after running liftover
gwas_stats_38 <- "./data/gwas/TRAF1_region_hg38.bed" |>
    read_tsv(col_names = c("chr", "start", "end", "rsid")) |>
    select(chr, pos = end, rsid) |>
    left_join(gwas_stats, by = c("rsid" = "SNP"))

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

# ref panel
Sys.setenv(VROOM_CONNECTION_SIZE = 500000L)

panel <- read_tsv("./data/chr9:119752571-121196432_JIApositions.vcf.gz", comment = "##")
plink <- read_delim("./data/chr9:119752571-121196432.ld", delim = " ", col_names = FALSE)
plink_filt <- plink[, sapply(plink, function(x) !all(is.na(x)))]
ldmat <- data.matrix(plink_filt)
rownames(ldmat) <- colnames(ldmat) <- panel$ID

freqx <- read_tsv("./data/chr9:119752571-121196432.frqx") |>
    select(CHR, ID = SNP, A1, A2, hom_A1 = 5, het = 6, hom_A2 = 7) |>
    mutate(ca1 = hom_A1*2 + het,
	   ca2 = hom_A2*2 + het,
	   fa1 = ca1 / (ca1 + ca2),
	   fa2 = ca2 / (ca1 + ca2)) |>
    select(CHR, ID, A1, A2, fa1, fa2) |>
    mutate(CHR = paste0("chr", CHR))

gwas_stats_38 |>
    select(chr, pos, major_allele, minor_allele) |>
    inner_join(select(panel, chr = `#CHROM`, pos = POS, ID, REF, ALT),
	       join_by(chr, pos)) |>
    mutate(REF_c = c("A" = "T", "T" = "A", "C" = "G", "G" = "C")[REF],
	   ALT_c = c("A" = "T", "T" = "A", "C" = "G", "G" = "C")[ALT]) |>
    mutate(m1 = major_allele == REF & minor_allele == ALT,
	   m2 = major_allele == ALT & minor_allele == REF,
	   m3 = major_allele == REF_c & minor_allele == ALT_c,
	   m4 = major_allele == ALT_c & minor_allele == REF_c) |>
    filter(ID == deviate$ID)

N <- stats_df |>
    filter(pos %in% panel$POS) |>
    mutate(n = n_case + n_ctrl) |>
    summarise(n = mean(n)) |>
    pull(n)

stats_df_filt <- stats_df |>
    filter(pos %in% panel$POS) |>
    mutate(pos = factor(pos, levels = panel$POS)) |>
    arrange(pos)

# run Susie
fit <- susie_rss(stats_df_filt$z, ldmat, n = N, coverage = 0.9)

png("./susie.png")
susie_plot(fit, y = 'PIP')
dev.off()

condz <- kriging_rss(stats_df_filt$z, ldmat, n = N, r_tol = 1e-04)
png("./susie_diagnostic.png")
condz$plot
dev.off()

cs_vars <- fit$sets$cs$L1

pip_df <- enframe(fit$pip, "ID", "pip") |>
    rowid_to_column() |>
    mutate(cs = ifelse(rowid %in% cs_vars, "CS1", "None")) |>
    extract(ID, "pos", "chr9_(\\d+)_.+", convert = TRUE, remove = FALSE) |>
    left_join(select(stats_df, pos, rsid, z)) |>
    mutate(z_exp = condz$conditional_dist$condmean)

write_tsv(pip_df, "./susie_pips.tsv")

deviate <- pip_df |> filter(abs(z - z_exp) == max(abs(z - z_exp)))

freqx |>
    filter(ID == deviate$ID)

panel |>
    filter(ID == deviate$ID)

