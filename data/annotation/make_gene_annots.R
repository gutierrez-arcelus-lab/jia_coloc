library(tidyverse)
library(readxl)

annotations <- 
    file.path("/lab-share/IM-Gutierrez-e2/Public/References/Annotations/hsapiens", 
              "gencode.v30.primary_assembly.annotation.gtf.gz") %>% 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
             col_names = c("chr", "feature", "start", "end", "strand", "info"))

gwas_stats <- read_excel("../gwas/TRAF1_region_summary_stats.xlsx")

gwas_stats_hg38 <- "../gwas/TRAF1_region_hg38.bed" %>%
    read_tsv(col_names = c("chr", "start", "end", "rsid")) %>%
    select(chr, pos = end, rsid) %>%
    left_join(gwas_stats, by = c("rsid" = "SNP"))

top_pos <- gwas_stats_hg38 %>%
    arrange(additive_pvalue) %>%
    slice(1) %>%
    pull(pos)

region <- annotations %>%
    filter(chr == "chr9", feature == "gene") %>%
    mutate(tss = case_when(strand == "+" ~ start, 
			   strand == "-" ~ end, 
			   TRUE ~ NA_integer_)) %>%
    filter(between(tss, top_pos - 2.5e5, top_pos + 2.5e5)) %>%
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+"),
	   gene_id = sub("^(ENSG\\d+)\\.\\d+((_PAR_Y)?)$", "\\1\\2", gene_id)) %>%
    select(chr, start, end, gene_id, gene_name)

gene_ids <- select(region, gene_id, gene_name)

write_tsv(gene_ids, "gene_annots.tsv")
