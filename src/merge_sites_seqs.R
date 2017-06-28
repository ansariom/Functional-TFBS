#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)

sites_file = args[1]
motif_seqs_file = args[2]
flank_seqs_file = args[3]
outfile = args[4]

d1 <- read.delim(sites_file, header = T, sep = "\t")

m <- read.table(motif_seqs_file, sep = "\t")
colnames(m) <- c("key", "motif_seq")

f <- read.table(flank_seqs_file, sep = "\t")
colnames(f) <- c("key", "flanking_seq")

d2 <- merge(m, f, by.x = "key", by.y="key")

d1$key <- paste(d1$chromosome, d1$tss_id, d1$feature, sep = "%")
d <- merge(d1, d2, by.x = "key", by.y="key")

d$key <- NULL

write.table(d, file = outfile, sep = "\t", quote = F, row.names = F)
