library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(ggsci)
library(ggtext)

myeloid_cells <- "monocyte|macrophage"
genes <- read_tsv("./data/coloc_inputs/gene_annots.tsv")

all_colocs <- read_tsv("./results/coloc_results_df_v2.tsv") |>
    filter(quant_method == "ge", study != "GTEx") |>
    drop_na(h4) |>
    select(study, tissue = tissue_label, condition = condition_label, 
	   group = qtl_group, gene_id, gene_name, pp4 = h4) |>
    filter(pp4 > 0.5) |>
    add_count(study, tissue, condition, gene_id) |>
    mutate(group = gsub("_", " ", group)) |>
    mutate(dataset = case_when(condition == "naive" & n == 1 ~ paste(study, tissue),
			       condition == "naive" & n > 1 ~ paste(study, group),
			       condition != "naive" ~ paste(study, tissue, condition))) |>
    select(dataset, gene_id, gene_name, pp4) |>
    mutate(dataset = gsub("_", " ", dataset),
	   mye = grepl(myeloid_cells, dataset),
	   cl = ifelse(mye == TRUE, ggsci::pal_npg()(10)[1], "grey20"),
	   gene_name = factor(gene_name, levels = sort(genes$gene_name)))

all_colocs$dataset[all_colocs$mye == TRUE] <- 
    paste0("<span style=\"color: ", 
	   all_colocs$cl[all_colocs$mye == TRUE] , 
	   "\">**", all_colocs$dataset[all_colocs$mye == TRUE], 
	   "**</span>")

all_colocs$dataset[all_colocs$mye == FALSE] <- 
    paste0("<span style=\"color: ", 
	   all_colocs$cl[all_colocs$mye == FALSE] , 
	   "\">", all_colocs$dataset[all_colocs$mye == FALSE], 
	   "</span>")
 

coloc_cols <- 
    c("white", "#a7b8c6", "white", "#fee17e", "#88c580", "white", 
      "mediumpurple1", "#f9c2c1", "#003f5c","#00a5ff", "tomato3") |>
    setNames(sort(genes$gene_name))

out <- all_colocs |>
    ggplot(aes(x = pp4, 
	       y = reorder(dataset, pp4, max))) +
    geom_point(aes(fill = gene_name),
	       size = 4, shape = 21, stroke = .3) +
    scale_fill_manual(values = coloc_cols, drop = FALSE) +
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

ggsave("./results/colocs_50.png", out, width = 8, height = 9, dpi = 600)


# Bayes factors per variant
#
#colocs_snp <- read_tsv("./results/coloc_results_df.tsv")
#
#colocs_snp |>
#    filter(study == "BLUEPRINT", tissue == "monocyte", gene_name == "TRAF1", feature == "ge") |>
#    select(study, tissue, condition, gene_name, snp, pp4 = SNP.PP.H4) |>
#    arrange(desc(pp4))
#    
#
