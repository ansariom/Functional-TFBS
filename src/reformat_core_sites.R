#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

intable = args[1]
orig_data <- read.delim(intable, sep = "\t", header = T)

orig_data <- orig_data[,c("chromosome", "strand", "transcript_id", "gene_id", "symbol", "full_name", "tss_mode_loc", "tss_mode_prob", "mode_read_count","peak_start", "peak_end", "shape", "peak_read_count", "model_W", "loglik_score", "pwm_name", "win_strand", "win_no", "rel_loc", "site_start", "site_end", "feature", "tss_id", "promoter_start",  "pwm_length")]

write.table(file=intable, orig_data, sep = "\t", col.names = T, row.names = F, quote = F)
