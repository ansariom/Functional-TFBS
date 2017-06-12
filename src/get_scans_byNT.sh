#!/bin/bash

# This script runs TFBS scanner for every individual NT in parallel
tfbs_scan_jar=/nfs0/BPP/Megraw_Lab/mitra/Projects/TFBSScanSuite/tfbs_scan.jar

input_fasta=$1
roe_fwd=$2
roe_rev=$3
pwm=$4
nucsAfterTSS=$5
ncpu=$6
outdir=$7
feature_file=$8

BGWin=250

if [ ! -d $outdir ];then
        mkdir -p $outdir
fi

# convert fasta format to oneline
#main_seq_fasta=$input_fasta.oneLine
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $input_fasta > $main_seq_fasta


seqs_outdir=$outdir/seqs
if [ ! -d $seqs_outdir ]; then
        mkdir -p $seqs_outdir
fi
seq_prefix=$seqs_outdir/seq.
nseqs=200

split -d -l $nseqs $input_fasta $seq_prefix

count=1
for seqfile in `ls $seqs_outdir/seq*`; do
	if [ $count -lt $ncpu ]; then
		outfile=$outdir/`basename $seqfile`.Rdat
		java -Xms100G -Xmx200G -jar $tfbs_scan_jar GenFeaturesByNT -N $nucsAfterTSS -B $BGWin $roe_fwd $roe_rev $seqfile $pwm $outfile &
		let "count=count+1"
        else
                echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
	fi
done
wait

echo "All finished!!!"

cat $seq_prefix*.Rdat | grep -v "_FWD_" > $outdir/tmp.Rdat
head -1 $seq_prefix"00.Rdat" | cat - $outdir/tmp.Rdat > $feature_file

rm -f $outdir/tmp.Rdat
rm -f -r $outdir

