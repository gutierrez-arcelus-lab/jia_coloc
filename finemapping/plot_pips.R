library(tidyverse)
library(ggrepel)

annot <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens/",
	      "gencode.v38.primary_assembly.annotation.gtf.gz") |>
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info"))

genes <- c("RP11-27I1.6", "C5", "CNTRL", "FBXW2", "MEGF9", "PHF19", "PSMD5", "TRAF1")

tx_support <- 
    annot |>
    filter(chr == "chr9", feature == "transcript") |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+")) |>
    filter(gene_name %in% genes) |>
    mutate(tsl = str_extract(info, "(?<=transcript_support_level\\s\")[0-9]"),
	   tsl = as.integer(tsl),
	   level = str_extract(info, "(?<= level )[0-9]"),
	   level = as.integer(level)) |>
    select(gene_id, gene_name, transcript_id, start, end, strand, tsl, level) |>
    mutate(len = map2_int(start, end, function(x, y) length(x:y))) |> 
    group_by(gene_id) |>
    filter(is.na(tsl) | tsl == min(tsl)) |>
    ungroup() |>
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id)) |>
    select(gene_id, gene_name, transcript_id, start, end, strand, len)

exon_annot <- annot |>
    filter(chr == "chr9", feature == "exon") |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+")) |>
    select(chr, start, end, strand, gene_id, gene_name, transcript_id) |>
    mutate(transcript_id = sub("\\.\\d+$", "", transcript_id)) |>
    inner_join(select(tx_support, transcript_id, len)) |>
    group_by(gene_id) |>
    filter(len == max(len)) |>
    ungroup() |>
    select(-len)

intron_annot <- exon_annot |>
    group_by(gene_id) |>
    mutate(start_gene = min(start), end_gene = max(end)) |>
    ungroup() |>
    arrange(chr, start_gene, start) |>
    group_by(gene_id, gene_name, transcript_id, start_gene, end_gene) |>
    reframe(data = bind_cols(tibble(start = c(unique(start_gene), end)),
			       tibble(end = c(start, unique(end_gene))))) |>
    ungroup() |>
    unnest(cols = c(data)) |>
    mutate(start_i = ifelse(start == end | start == start_gene, start, start + 1L),
	   end_i = ifelse(start == end | end == end_gene, end, end - 1L)) |>
    filter(start != end) |>
    select(gene_id, gene_name, transcript_id, start = start_i, end = end_i)

annot_out <- 
    bind_rows("exon" = select(exon_annot, gene_id, gene_name, transcript_id, start, end),
              "intron" = intron_annot,
              .id = "feature") |>
    arrange(gene_name, gene_id, transcript_id, start) |>
    group_by(transcript_id) |>
    mutate(i = row_number(),
	   col = case_when(feature == "intron" & (i == min(i) | i == max(i)) ~ NA_character_,
			   TRUE ~ feature)) |>
    ungroup() |>
    filter(!is.na(col)) |>
    arrange(start, gene_id, gene_name, transcript_id) |>
    mutate(ytx = as.numeric(fct_inorder(transcript_id))) |>
    mutate(ygene = ifelse(gene_id == "ENSG00000270917.1", 0.07, 0.065))
    

arrows_df <- annot_out |>
    group_by(transcript_id, ytx) |>
    summarise(left = min(start), right = max(end)) |>
    mutate(x = map2(left, right, function(x, y) x:y)) |>
    select(transcript_id, ytx, x) |>
    unnest(cols = x) |>
    mutate(i = cut_interval(x, length = 1e4)) |>
    group_by(transcript_id, ytx, i) |>
    filter(row_number() == 1) |>
    left_join(select(tx_support, transcript_id, strand)) |>
    group_by(transcript_id) |>
    filter(n() == 1 | 
	   (n() > 1 & strand == "+" & row_number() != first(row_number())) |
	   (n() > 1 & strand == "-" & row_number() != last(row_number()))) |>
    ungroup() |>
    left_join(distinct(exon_annot, transcript_id, strand)) |>
    mutate(x2 = ifelse(strand == "+", x + 50, x - 50)) |>
    left_join(select(tx_support, gene_name, transcript_id)) |>
    mutate(ygene = ifelse(gene_name == "RP11-27I1.6", 0.07, 0.065))

gene_lab <- annot_out |>
    group_by(gene_id, transcript_id) |>
    slice(which.min(i)) |>
    ungroup() |>
    select(gene_name_19 = gene_name, transcript_id, ytx, pos = start) |>
    left_join(select(tx_support, gene_name, transcript_id)) |>
    mutate(ygene = ifelse(gene_name == "RP11-27I1.6", 0.07, 0.065))

pips <- read_tsv("./susie_pips.tsv")


out <- 
    ggplot(pips |> filter(pos > 120400000, cs == "None"), 
	   aes(pos, pip)) +
    geom_point() +
    geom_point(data = filter(pips, cs != "None"), color = "#de425b") +
    geom_text_repel(data = filter(pips, rsid %in% c("rs7039505", "rs7034653")),
		    aes(pos, pip, label = rsid),
		    segment.size = .1, min.segment.length = 0, 
		    size = 2) +
    geom_segment(data = arrows_df,
		 aes(x = x, xend = x2, y = ygene, yend = ygene),
		 linewidth = .2,
		 arrow = arrow(length = unit(0.1, "cm")),
		 color = "midnightblue") +
    geom_segment(data = filter(annot_out, feature == "intron"),
		 aes(x = end, xend = start, y = ygene, yend = ygene),
		 linewidth = .5, 
		 color = "midnightblue") +
    geom_segment(data = filter(annot_out, feature == "exon"),
		 aes(x = start, xend = end, y = ygene, yend = ygene),
		 linewidth = 4, 
		 color = "midnightblue") +
    geom_text(data = gene_lab, 
	      aes(x = pos, y = ygene, label = gene_name),
	      fontface = "italic", 
	      size = 2,
	      hjust = 0,
	      nudge_y = c(rep(0.0035, 7), 0.0035),
	      color = "grey35") +
    scale_x_continuous(labels = function(x) x/1e6L) +
    theme(panel.background = element_rect(fill = "grey96")) +
    labs(x = "Position in chr9 (Mb, GRCh38)", 
	 y = "Posterior Inclusion Probability")

ggsave("./pips.png", out, height = 3, width = 6)

