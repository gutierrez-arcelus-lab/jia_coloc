library(tidyverse)
library(RColorBrewer)

coloc_h <- read_tsv("./coloc_results_summary.tsv")

coloc_h_plotdf <- coloc_h %>%
    filter(study != "GTEx") %>%
    mutate(dataset = paste0(study, " (", qtl_group, ")")) %>%
    select(dataset, gene_name, feature, h4) %>%
    group_by(dataset, feature) %>%
    filter(any(h4 > 0.8)) %>%
    ungroup() %>%
    group_by(dataset, feature, gene_name) %>%
    filter(h4 == max(h4)) %>%
    ungroup() %>%
    mutate(feature = recode(feature, "microarray" = "array", "ge" = "gene expression"),
	   feature = factor(feature, levels = c("gene expression", "array", "exon", "tx", "txrev")))

coloc_cols <- rev(c(ggsci::pal_npg()(10), "black"))
coloc_cols[4] <- "goldenrod3"

out <- ggplot(data = coloc_h_plotdf, 
       aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    facet_grid(feature~., scale = "free", space = "free_y") +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  plot.caption = element_text(hjust = 0)) +
    labs(x = "Probability of shared causal signal", 
	 y = NULL, 
	 fill = "Gene",
	 title = "eQTL catalogue datasets with any evidence of\ncolocalization for the TRAF1-C5 locus in JIA",
	 subtitle = "PP4 > 0.8 was only observed for eQTLs")

ggsave("./pp4.png", out, height = 14, width = 10)

coloc_snp <- read_tsv("./coloc_results_df.tsv")

snp_pp4_df <- coloc_h_plotdf %>% 
    filter(gene_name == "TRAF1", h4 > 0.75) %>%
    extract(dataset, c("study", "qtl_group"), "(.+) \\((.+)\\)", remove = FALSE) %>%
    select(dataset:feature) %>%
    inner_join(coloc_snp) %>%
    group_by(snp) %>%
    filter(any(SNP.PP.H4 > .05)) %>%
    ungroup() %>%
    select(dataset, study, qtl_group, snp, snp_pp4 = SNP.PP.H4)

arrange(snp_pp4_df, desc(snp_pp4))

out2 <- 
    ggplot(snp_pp4, df, 
	   aes(x = reorder(snp, snp_pp4, mean), y = snp_pp4)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = NULL)

ggsave("./pp4_perSNP.png", out2, width = 6, height = 3)


