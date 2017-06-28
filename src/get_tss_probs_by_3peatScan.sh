#!/bin/bash

model_file=$1
roe_fwd_tbl=$2
roe_rev_tbl=$3
peat_pwms=$4

promoter_hw=$5  # The half width of input sequence fasta, centered at 0 (TSS) (such as -5000 to 5000 so set it to 5000)
scan_hw=$6	# How many nts we want to scan for 3PEAT? Give the half width like (-10 to 10 so set it to 10)
scripts_dir=software

l1logreg_path=/nfs0/BPP/Megraw_Lab/mitra/software/3PEAT_Model/l1_logreg-0.8.2-i686-pc-linux-gnu

### Paralle processing for many FASTA seqs
###########################################

nseqs=100  # number of seqs in each inputfile for parallel execution
let "nlocs=$nseqs/2"

outdir=$7
input_fasta=$8
input_loc_bed=$9
outfile=${10}

if [ ! -d $outdir ];then
        mkdir -p $outdir
fi

seqs_outdir=$outdir/Sequences
if [ ! -d $seqs_outdir ]; then
        mkdir -p $seqs_outdir
fi
seq_prefix=$seqs_outdir/seq.
loc_prefix=$seqs_outdir/loc.

split -d -l $nseqs $input_fasta $seq_prefix
split -d -l $nlocs $input_loc_bed $loc_prefix


count=1
ncpu=30
for seqfile in `ls $seqs_outdir/seq*`; do
        if [ $count -lt $ncpu ]; then
                filename=$(basename "$seqfile")
                i="${filename##*.}"
                bedfile=$loc_prefix"$i"
                echo $bedfile
                echo $seqfile
		echo software/3PEAT_Scan.sh $model_file $seqfile $bedfile $peat_pwms $roe_fwd_tbl $roe_rev_tbl $scan_hw $promoter_hw $outdir software $l1logreg_path
		software/3PEAT_Scan.sh $model_file $seqfile $bedfile $peat_pwms $roe_fwd_tbl $roe_rev_tbl $scan_hw $promoter_hw $outdir software $l1logreg_path
                let "count=count+1"
        else
                echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
        fi
done
wait

echo "All finished!!!"

rm -f $outfile

for wig_file in  $outdir/scans/seq.*/*.wig; do
	peak_file=`basename $wig_file`
        peak_name=${peak_file%.*}
	software/extract_prob_at_tss.pl $wig_file $peak_name >> $outfile
done

rm -f -r $seqs_outdir





