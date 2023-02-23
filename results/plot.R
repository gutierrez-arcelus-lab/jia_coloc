library(tidyverse)
library(tidytext)
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
    mutate(study = sub("_", " ", study),
	   qtl_group = gsub("_", " ", qtl_group),
	   dataset = paste0(study, " (", qtl_group, ")"),
	   feature = recode(feature, "ge" = "RNA-seq", "microarray" = "Microarray"),
	   feature = factor(feature, levels = c("RNA-seq", "Microarray"))) %>%
    select(dataset, gene_name, feature, h4) %>%
    group_by(dataset, feature, gene_name) %>%
    filter(h4 == max(h4)) %>%
    ungroup()

coloc_cols <- c("black", ggsci::pal_npg()(10))
coloc_cols <- rev(coloc_cols[c(1, 2, 3, 6, 7, 8, 11)])
names(coloc_cols) <- sort(unique(coloc_mono$gene_name))

fig_mono <- coloc_mono %>%
    ggplot(aes(x = h4, y = dataset, fill = gene_name)) +
    geom_vline(xintercept = 0.8, linetype = 2, size = .25) +
    geom_point(size = 3, shape = 21, stroke = .2) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    facet_grid(feature~., scale = "free_y", space = "free", switch = "y") +
    theme(text = element_text(size = 12),
	  panel.background = element_rect(fill = "grey96"),
	  strip.placement = "outside") +
    labs(x = "Probability of shared signal", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./fig_mono.png", fig_mono, height = 5, width = 6, dpi = 600)
ggsave("./fig_mono.tiff", fig_mono, height = 5, width = 6, dpi = 600)



### Supplement figure

immune_cells <- c("monocyte|macrophage")

read_tsv("./coloc_results_df_v2.tsv") |>
    filter(quant_method == "ge", study != "GTEx") |>
    distinct(study, tissue_label, condition_label, qtl_group)


all_colocs <- read_tsv("./coloc_results_df_v2.tsv") |>
    filter(quant_method == "ge", study != "GTEx") |>
    drop_na(h4) |>
    select(study, tissue = tissue_label, condition = condition_label, 
	   group = qtl_group, gene_id, gene_name, pp4 = h4) |>
    group_by(study, tissue, condition, group, gene_id, gene_name) |>
    slice(which.max(pp4)) |>
    ungroup() |>
    filter(pp4 > 0.5) |>
    add_count(study, tissue, condition, gene_id) |>
    mutate(group = gsub("_", " ", group)) |>
    mutate(dataset = case_when(condition == "naive" & n == 1 ~ paste(study, tissue),
			       condition == "naive" & n > 1 ~ paste(study, tissue, " - ", group),
			       condition != "naive" ~ paste0(study, " ", tissue, " ", condition))) |>
    select(dataset, gene_id, gene_name, pp4) |>
    mutate(dataset = gsub("_", " ", dataset)) |>
    mutate(immune = grepl(immune_cells, dataset),
	   cl = ifelse(immune == TRUE, ggsci::pal_npg()(10)[1], "grey20"))

all_colocs$dataset[all_colocs$immune == TRUE] <- 
    paste0("<span style=\"color: ", 
	   all_colocs$cl[all_colocs$immune == TRUE] , 
	   "\">**", all_colocs$dataset[all_colocs$immune == TRUE], 
	   "**</span>")

all_colocs$dataset[all_colocs$immune == FALSE] <- 
    paste0("<span style=\"color: ", 
	   all_colocs$cl[all_colocs$immune == FALSE] , 
	   "\">", all_colocs$dataset[all_colocs$immune == FALSE], 
	   "</span>")
 
sort(unique(all_colocs$gene_name))

coloc_cols <- c("AL161911.1" = "lightskyblue", 
		"C5" = "#f9d14a",    
		"CNTRL" = "lightgreen",
		"FBXW2" = "grey55",
		"MEGF9" = "lightpink",
		"PHF19" = "royalblue3",
		"PSMD5" = "white", 
		"TRAF1" = "tomato3")

out2 <- all_colocs |>
    ggplot(aes(x = pp4, 
	       y = reorder(dataset, pp4, max))) +
    geom_point(aes(fill = gene_name),
	       size = 4, shape = 21, stroke = .3) +
    scale_fill_manual(values = coloc_cols) +
    scale_x_continuous(breaks = seq(0.5, 1, by = .1)) +
    theme(text = element_text(size = 12),
	  axis.text.x = element_text(size = 12), 
	  axis.title.x = element_text(size = 12), 
	  axis.text.y = ggtext::element_markdown(size = 11),
	  panel.grid.minor.x = element_blank(),
	  panel.background = element_rect(fill = "grey95")) +
    labs(x = "Posterior probability of colocalization", 
	 y = NULL, 
	 fill = "Gene")

ggsave("./colocs_50.png", out2, width = 7, height = 9, dpi = 600)

greens <- grep("green", colors(), value=T)

png("test.png")
plot(1:40, 1:40, col = greens, pch = 19, cex = 2)
dev.off()


