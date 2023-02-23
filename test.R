library(tidyverse)
library(readxl)
library(ggrepel)
library(coloc)

stats <- read_excel("./data/gwas/TRAF1_region_summary_stats.xlsx")
snps <- c("rs7039505", "rs7034653")

## hg19
gwas19 <- "./data/gwas/TRAF1_region_hg19.bed" %>%
    read_tsv(col_names = c("chr", "start", "end", "snp")) %>%
    select(pos = end, snp) %>%
    left_join(stats, by = c("snp" = "SNP")) %>%
    select(snp, pos, everything())

hit_pos19 <- gwas19 %>% 
    arrange(additive_pvalue) %>% 
    slice(1) %>%
    pull(pos)

gwas19_trim <- filter(gwas19, between(pos, hit_pos19 - 5e5, hit_pos19 + 5e5))

## compare with open GWAS harmonized data
gwas_df <- gwas19_trim %>%
    select(snp, position = pos, n_case, n_ctrl, or = additive_odds_ratio, pvalue = additive_pvalue) %>%
    mutate(beta = log(1/or),
	   varbeta = (beta^2) / qchisq(pvalue, df = 1, lower = FALSE),
	   type = "cc",
	   s = n_case/n_ctrl) %>%
    select(beta, varbeta, snp, position, type, s, pvalue)

stats_harmonized <- "./data/gwas/ebi-a-GCST005528.vcf.gz" %>%
    read_tsv(comment = "##")

stats_harmonized_trim <- stats_harmonized %>%
    filter(`#CHROM` == 9, between(POS, hit_pos19 - 5e5, hit_pos19 + 5e5)) %>%
    separate(`EBI-a-GCST005528`, c("effsize", "se", "log_p", "id"), sep = ":", convert = TRUE) %>%
    select(position = POS, snp = ID, effsize, se, log_p)

gwas_catalog <- "./data/gwas/23603761-GCST005528-EFO_1001999-Build37.f.tsv.gz" %>%
    read_tsv() %>%
    filter(chromosome == 9) %>%
    select(position = base_pair_location, snp = variant_id, log_p = p_value, beta, se = standard_error) %>%
    mutate(log_p = -log10(log_p))

gwas_catalog_trim <- gwas_catalog %>%
    filter(between(position, hit_pos19 - 5e5, hit_pos19 + 5e5))

gwas_plot_df <- 
    bind_rows("GWAS catalog" = gwas_catalog_trim,
	      "New Summary stats" = select(gwas_df, snp, position, log_p = pvalue) %>%
		  mutate(log_p = -log10(log_p)),
	      .id = "source")

plot1 <- ggplot(gwas_plot_df, aes(x = position, y = log_p)) +
    geom_hline(yintercept = -log10(1e-6), color = "grey") +
    geom_hline(yintercept = -log10(5e-8), color = "black") +
    geom_point(color = "grey60") +
    geom_point(data = filter(gwas_plot_df, snp %in% snps)) +
    geom_text_repel(data = filter(gwas_plot_df, snp %in% snps), 
		    aes(label = snp),
		    min.segment.length = 0) +
    scale_x_continuous(labels = function(x) x/1e6) +
    facet_wrap(~source, ncol = 1) +
    theme(panel.background = element_rect(fill = "grey98")) +
    labs(x = "Position in chr9 (Mb)")

plot2 <- gwas_df %>%
    inner_join(gwas_catalog_trim, by = c("snp", "position")) %>%
    ggplot(aes(beta.x, beta.y)) +
    geom_point() +
    labs(x = "Beta (new stats)", y = "Beta (catalog)")

plot3 <- gwas_df %>%
    inner_join(gwas_catalog_trim, by = c("snp", "position")) %>%
    ggplot(aes(abs(beta.x), abs(beta.y))) +
    geom_point() +
    labs(x = "Abs beta (new stats)", y = "Abs beta (catalog)")

plot4 <- gwas_df %>%
    inner_join(gwas_catalog_trim, by = c("snp", "position")) %>%
    ggplot(aes(varbeta, se^2)) +
    geom_point() +
    labs(x = "Var beta (new stats)", 
	 y = "Var beta (catalog)")

plot5 <- gwas_df %>%
    inner_join(gwas_catalog_trim, by = c("snp", "position")) %>%
    ggplot(aes(-log10(pvalue), log_p)) +
    geom_point() +
    labs(x = "-log10 (pvalue) (new stats) ", 
	 y = "-log10 (pvalue) (catalog)")

out <- 
    cowplot::plot_grid(plot1,
		       NULL,
		       cowplot::plot_grid(plot2, plot3, plot4, plot5, ncol = 2),
		       ncol = 1,
		       rel_heights = c(1, .05, 1)) + 
    theme(panel.background = element_rect(fill = "white", color = "white"))

ggsave("./jia_gwas.png", out, height = 8)


# check allele order

gwas_catalog_gt <- "./data/gwas/23603761-GCST005528-EFO_1001999-Build37.f.tsv.gz" %>%
    read_tsv() %>%
    filter(chromosome == 9) %>%
    select(pos = base_pair_location, rsid = variant_id, other_allele, effect_allele)

gwas_new_gt <- gwas19_trim %>%
    select(pos, rsid = snp, minor_allele, major_allele, new_maf = maf_ctrl)

gwas_harm_gt <- stats_harmonized %>%
    select(pos = POS, rsid = ID, REF, ALT)

left_join(gwas_new_gt, gwas_catalog_gt) %>%
    left_join(gwas_harm_gt) %>%
    select(pos, rsid, new_maf, new_major = major_allele, new_minor = minor_allele,  
	   cat_other = other_allele, cat_effect = effect_allele,
	   hd_ref = REF, hd_alt = ALT) %>%
    filter(new_minor != cat_effect)
