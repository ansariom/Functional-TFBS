#!/bin/bash

sites_table=$1
genome=$2
motifs_fasta=$3
flanking_len=$4

#cat $sites_table | awk -F "\t" '(NR > 1){chr=$1; strand=$2; start=$20-1; end=$21-1; name=chr"%"$23"%"$22; print chr"\t"start"\t"end"\t"name"\t.\t"strand}' > $sites_table.bed

# correction on site chromosomal location: 
# The relative locations in raw_loglik scores are computed based on rev-comp of -strands for antisense genes. So 
# Since we get getfasta from oriinal genome coordinates we need to correct the start/end locarions for sites in - strand

motifs_bed=$sites_table.bed

cat $sites_table | awk -F "\t" '(NR > 1){chr=$1; strand=$2; mode = $7; site_start=$20; rel_loc = $19; tss_id = $23; pwm_len = $25; feature = $22; if (strand == "-") {end = mode - (rel_loc); start = end - pwm_len  }  else {start = site_start - 1; end = start + pwm_len }; name=chr"%"tss_id"%"feature; print chr"\t"start"\t"end"\t"name"\t.\t"strand}' > $motifs_bed

# get the motif (PWM binding) sequences according to ROE_win_locs relative to tss mode.
motifs_out=$motifs_bed.fa
bedtools getfasta -fi $genome -bed $motifs_bed -fo $motifs_out -name -tab -s


flanking_bed=$sites_table.flanking.bed
flanking_out=$sites_table.flanking.bed.fa
cat $motifs_bed | awk -v f=$flanking_len '{start = $2 - f; end = $3 + f; print $1"\t"start"\t"end"\t"$4"\t"$5"\t"$6 }' > $flanking_bed
# get the sequnces in vicinity of motis_seq within -flanking_len to +flanking_len
bedtools getfasta -fi $genome -bed $flanking_bed -fo $flanking_out -name -tab -s

# merge the two sequences file into final output
software/merge_sites_seqs.R $sites_table $motifs_out $flanking_out $motifs_fasta
