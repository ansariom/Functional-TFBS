#!/usr/bin/Rscript

process_frame <- function(d, fileFrame) {
        d$tss_id <- rownames(fileFrame)
        dlong <- gather(d, feature, loglik_score, -tss_id)
        dlong <- subset(dlong, dlong$loglik_score > 0)

        # get PWM name as separate column
        dlong$pwm_name <- lapply(dlong$feature, function(x) { gsub("_FWD_\\d+_\\D*\\d+", "", x) })
        dlong$pwm_name <- lapply(dlong$pwm_name, function(x) { gsub("_REV_\\d+_\\D*\\d+", "", x) })

        # get window number as separate column
        dlong$win_no <- lapply(dlong$feature, function(x) { gsub("[[:graph:]]*_FWD_(\\d+)_\\D*\\d+", "\\1", x) })
        dlong$win_no <- lapply(dlong$win_no, function(x) { gsub("[[:graph:]]*_REV_(\\d+)_\\D*\\d+", "\\1", x) })

        # get the location relative to TSS
        dlong$rel_loc <- lapply(dlong$feature, function(x) { gsub("[[:graph:]]*_FWD_\\d+_(\\D*\\d+)", "\\1", x) })
        dlong$rel_loc <- lapply(dlong$rel_loc, function(x) { gsub("[[:graph:]]*_REV_\\d+_(\\D*\\d+)", "\\1", x) })

        dlong$win_strand <- lapply(dlong$feature, function(x) { gsub("[[:graph:]]*_(FWD)_\\d+_\\D*\\d+", "\\1", x) })
        dlong$win_strand <- lapply(dlong$win_strand, function(x) { gsub("[[:graph:]]*_(REV)_\\d+_\\D*\\d+", "\\1", x) })

        # get the tss info
        dlong$sepcol <- dlong$tss_id
        dlong <- separate(dlong, sepcol, sep = "_", into = c("transcript_id", "chromosome", "promoter_start", "t"))
        dlong$t <- NULL

        dlong$sepcol <- dlong$transcript_id
        dlong <- separate(dlong, sepcol, sep = "\\.", into = c("gene_id", "index"))
        dlong$index <- NULL

        dlong$tss_mode_loc <- as.numeric(dlong$promoter_start) + nucs_upstream
        dlong$site_start <- dlong$tss_mode_loc + as.numeric(dlong$rel_loc)

        dlong$pwm_name <- unlist(dlong$pwm_name)
        dlong$win_no <- unlist(dlong$win_no)
        dlong$rel_loc <- unlist(dlong$rel_loc)
	dlong$win_strand <- unlist(dlong$win_strand)

        # Merge with PWM length info to get site coordinates
        dlong <- merge(dlong, pwms_length, by.x = "pwm_name", by.y = "pwm")
        dlong$site_end <- dlong$site_start + dlong$pwm_length

        # 2- merge with other infrmation
        dlong <- merge(dlong, genes, by.x = "gene_id", by.y = "locus_name", all.x = T)
        dlong <- merge(dlong, peaks, by.x = c("tss_mode_loc", "gene_id"), by.y = c("ModeLocation", "GeneName"))
        drops <- c("Chromosome", "ModeLocation", "TranscriptLocation", "TranscriptID", "GeneName","GeneType","PeakID")
        dlong <- dlong[, !(names(dlong) %in% drops)]

        # merge with tss probs at TSS mode location
        dlong <- merge(dlong, tss_probs, by.x = "tss_id", by.y = "key", all.x = T)

        # Merge with model coeffs and compute the total score
        dlong$pwm_win <- paste(dlong$pwm_name, dlong$win_strand, dlong$win_no, sep = "_")
        dlong <- merge(dlong, model_w, by.x = c("pwm_win"), by.y = c("pwm_win"))
        # compute totalScore
        dlong$total_score <- dlong$model_W*dlong$loglik_score
        dlong$tss_name <- NULL
	dlong$pwm_win <- NULL

        colnames(dlong) <- c("tss_id", "tss_mode_loc", "gene_id", "pwm_name", "feature", "loglik_score", "win_no", "rel_loc", "win_strand" ,"transcript_id", "chromosome", "promoter_start", "site_start", "pwm_length", "site_end", "symbol", "full_name", "strand", "peak_start", "peak_end", "peak_read_count", "mode_read_count", "shape", "tss_mode_prob", "model_W", "total_score")
	return(dlong)
}

library(tidyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1] # input rdat file containing loglik scores for each nt. Header contains feature names (pwm_(FWD/REV)_WIN_(LOC)
gene_info_file <- args[2]
peaks_file <- args[3]
pwm_length_file <- args[4]
tss_prob_file <- args[5]
model_weights <- args[6]
nucs_upstream <- as.numeric(args[7])
outfile <- args[8]

fileFrame = fread(infile, header=TRUE)
peaks <- read.csv(peaks_file, header = T)
genes <- read.delim(gene_info_file, header = T, sep = "\t")
pwms_length <- read.table(pwm_length_file, header = F)
tss_probs <- read.table(tss_prob_file, header = F)
model_w <- read.table(model_weights, header = F)

colnames(pwms_length) <- c("pwm", "pwm_length")

colnames(tss_probs) <- c("tss_name", "tss_mode_prob")

colnames(model_w) <- c("pwm_win" , "model_W")
# sort based on weights pick top 50
model_w <- model_w[order(-abs(model_w$model_W)),]
model_w <- model_w[1:50,]

# add key to tss_probs
tss_probs$key <- paste(tss_probs$tss_name, "_0", sep = "")


# 1 - Reshape the table
# The reshaping of large columns is taking too long so we broke down the process into smaller parts to be faster in process
rownames(fileFrame) <- fileFrame$V1
fileFrame$V1 <- NULL
fileFrame$tss_id <- rownames(fileFrame)

# break down the data into subset of columns (columns are amny ~123000)
ncols <- ncol(fileFrame)
p <- 10000
parts <- floor(ncols/p)
start = 1

all <- data.frame(tss_id=c(), tss_mode_loc=c(), gene_id=c(), pwm_name=c(), feature=c(), loglik_score=c(), win_no=c(), rel_loc=c(), win_strand=c() ,transcript_id=c(), chromosome=c(), promoter_start=c(), site_start=c(), pwm_length=c(), site_end=c(), symbol=c(), full_name=c(), strand=c(), peak_start=c(), peak_end=c(), peak_read_count=c(), mode_read_count=c(), shape=c(), tss_mode_prob=c(), model_W=c(), total_score=c())
#all <- data.frame(tss_id=c(), feature=c(), loglik_score=c(), pwm_name=c())
for(i in 1:parts) {
	end <- i * p
	d <- fileFrame[, start:end]
	all <- rbind(all, process_frame(d, fileFrame))
	start <- end + 1
}

# run for remaining items
d <- fileFrame[, start:ncols]
all <- rbind(all, process_frame(d, fileFrame))

write.table(file=outfile, all, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t");
