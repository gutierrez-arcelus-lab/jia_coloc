library(tidyverse)
library(readxl)
library(ggrepel)
library(coloc)

stats <- read_excel("./data/gwas/TRAF1_region_summary_stats.xlsx")
snps <- c("rs7039505", "rs7034653")

gwas38 <- "./data/gwas/TRAF1_region_hg38.bed" %>%
    read_tsv(col_names = c("chr", "start", "end", "snp")) %>%
    select(pos = end, snp) %>%
    left_join(stats, by = c("snp" = "SNP")) %>%
    select(snp, pos, everything())

hit_pos <- gwas38 %>% 
    arrange(additive_pvalue) %>% 
    slice(1) %>%
    pull(pos)

gwas38_trim <- filter(gwas38, between(pos, hit_pos - 5e5, hit_pos + 5e5))

ggplot(gwas38_trim, aes(x = pos, y = -log10(best_pvalue))) +
    geom_point(color = "grey70") +
    geom_point(data = filter(gwas38_trim, snp %in% snps)) +
    geom_text_repel(data = filter(gwas38_trim, snp %in% snps), 
		    aes(label = snp),
		    min.segment.length = 0) +
    theme_bw()

minimum_df <- gwas38_trim %>%
    select(snp, position = pos, n_case, n_ctrl, or = additive_odds_ratio, pvalue = additive_pvalue) %>%
    mutate(beta = ifelse(or >= 1, log(or), log(1/or)),
	   varbeta = (beta^2) / qchisq(pvalue, df = 1, lower = FALSE),
	   type = "cc",
	   s = n_case/n_ctrl) %>%
    select(beta, varbeta, snp, position, type, s)

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

minimum_df_hg19 <- gwas19_trim %>%
    select(snp, position = pos, n_case, n_ctrl, or = additive_odds_ratio, pvalue = additive_pvalue) %>%
    mutate(beta = ifelse(or >= 1, log(or), log(1/or)),
	   varbeta = (beta^2) / qchisq(pvalue, df = 1, lower = FALSE),
	   type = "cc",
	   s = n_case/n_ctrl) %>%
    select(beta, varbeta, snp, position, type, s)


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


gwas_ <- select(gwas19_trim, rsid = snp, ref_allele, major_allele, minor_allele, maf_ctrl)

stats_harmonized %>%
    select(rsid = ID, REF, ALT) %>%
    inner_join(gwas_) %>%
    arrange(desc(maf_ctrl))





gwas_plot_df <- 
    bind_rows(
	      "OpenGWAS" = select(stats_harmonized_trim, snp, position, log_p),
	      "Qiang" = select(gwas_df, snp, position, log_p = pvalue) %>%
		  mutate(log_p = -log10(log_p)),
	      .id = "source")

plot1 <- ggplot(gwas_plot_df, aes(x = position, y = log_p)) +
    geom_point(color = "grey60") +
    geom_point(data = filter(gwas_plot_df, snp %in% snps)) +
    geom_text_repel(data = filter(gwas_plot_df, snp %in% snps), 
		    aes(label = snp),
		    min.segment.length = 0) +
    facet_wrap(~source, scale = "free_y", ncol = 1) +
    theme(panel.background = element_rect(fill = "grey98"))

plot2 <- gwas_df %>%
    inner_join(stats_harmonized_trim, by = c("snp", "position")) %>%
    ggplot(aes(beta, effsize)) +
    geom_point() +
    labs(x = "Beta (from log(OR))", y = "Beta (Open GWAS)")

plot3 <- gwas_df %>%
    inner_join(stats_harmonized_trim, by = c("snp", "position")) %>%
    ggplot(aes(varbeta, se^2)) +
    geom_point() +
    labs(x = "Var (beta) (from log(OR)/p-value quantile)", 
	 y = "Var (Beta) (Open GWAS)")

plot4 <- gwas_df %>%
    inner_join(stats_harmonized_trim, by = c("snp", "position")) %>%
    ggplot(aes(-log10(pvalue), log_p)) +
    geom_point() +
    labs(x = "-log10 (pvalue) ", 
	 y = "p-value (Open GWAS)")

gwas_df %>%
    inner_join(stats_harmonized_trim, by = c("snp", "position")) %>%
    filter(se^2 > 0 & varbeta == 0) %>%
    arrange(desc(se))

out <- 
    cowplot::plot_grid(plot1, 
		       cowplot::plot_grid(plot2, plot3, plot4, ncol = 2),
		       ncol = 1) + 
    theme(panel.background = element_rect(fill = "white", color = "white"))
    

ggsave("./jia_gwas.png", out, height = 8)




