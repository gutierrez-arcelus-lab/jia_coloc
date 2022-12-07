library(tidyverse)
library(RColorBrewer)
library(ggsci)

coloc_h <- read_tsv("./coloc_results_summary.tsv")

coloc_rnaseq <- coloc_h %>%
    filter(study != "GTEx_V8", feature != "microarray") %>%
    mutate(dataset = paste0(study, " (", qtl_group, ")")) %>%
    select(dataset, gene_name, feature, h4) %>%
    group_by(dataset) %>%
    filter(any(h4 > 0.8)) %>%
    ungroup() %>%
    group_by(dataset, feature, gene_name) %>%
    filter(h4 == max(h4)) %>%
    ungroup() %>%
    mutate(feature = recode(feature, "ge" = "gene expression"),
	   feature = factor(feature, levels = c("gene expression", "exon", "tx", "txrev")))

coloc_array <- coloc_h %>%
    filter(feature == "microarray") %>%
    mutate(dataset = paste0(study, " (", qtl_group, ")")) %>%
    select(dataset, gene_name, feature, h4) %>%
    group_by(dataset) %>%
    filter(any(h4 > 0.8)) %>%
    ungroup() %>%
    group_by(dataset, feature, gene_name) %>%
    filter(h4 == max(h4)) %>%
    ungroup() %>%
    mutate(feature = recode(feature, "microarray" = "array"))

coloc_cols <- rev(c(ggsci::pal_npg()(10), "black"))
coloc_cols[4] <- "goldenrod3"
names(coloc_cols) <- sort(unique(coloc_h$gene_name))

out <- ggplot(data = coloc_rnaseq, 
       aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid(.~feature) +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./pp4.png", out, height = 10, width = 10)

out_array <- ggplot(data = coloc_array, 
       aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid(.~feature) +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  panel.border = element_rect(color = "black", fill = NA, size = 1),
	  legend.position = "none") +
    labs(x = NULL, 
	 y = NULL)

ggsave("./pp4_array.png", out_array, height = 2, width = 3)



out_mono <- coloc_rnaseq %>%
    filter(grepl("(macrophage)|(monocyte)", dataset, ignore.case = TRUE)) %>%
    ggplot(data = ., 
	   aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid(.~feature) +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./pp4_monocytes.png", out_mono, width = 10, height = 5)


out_mono_array <- coloc_array %>%
    filter(grepl("(macrophage)|(monocyte)", dataset, ignore.case = TRUE)) %>%
    ggplot(data = ., 
	   aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, color = "grey50") +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_grid(.~feature) +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  panel.border = element_rect(color = "black", fill = NA, size = 1),
	  legend.position = "none") +
    labs(x = NULL, 
	 y = NULL)

ggsave("./pp4_monocytes_array.png", out_mono_array, width = 3.3, height = 2)


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

# Monocytes for paper

coloc_mono <- coloc_h %>%
    filter(grepl("(monocyte|macrophage)", qtl_group, ignore.case = TRUE)) %>%
    filter(feature %in% c("ge", "microarray")) %>%
    mutate(dataset = paste0(study, " (", qtl_group, ")"),
	   dataset = ifelse(feature == "microarray", paste(dataset, "*"), dataset))%>%
    select(dataset, gene_name, feature, h4) %>%
    group_by(dataset, feature, gene_name) %>%
    filter(h4 == max(h4)) %>%
    ungroup()

coloc_cols <- rev(c("black", ggsci::pal_npg()(10)))
coloc_cols[3] <- "goldenrod3"
names(coloc_cols) <- sort(unique(coloc_h$gene_name))

fig_mono <- coloc_mono %>%
    ggplot(aes(x = h4, y = dataset, fill = gene_name)) +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96")) +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./fig_mono.png", fig_mono, height = 5, width = 6, dpi = 600)
