#!/bin/bash

# input: features_byNT.Rdat, tss_probs_at_mode, gene_aliases

# output a table with the following header

features_rdat=$1	#features_byNT.Rdat
outdir=$2
tss_probs=$3
pwms_file=$4
genes_aliases_file=$5 # gene alises
peaks_file=$6 # peaks.csv in david's peak caller format
model_weights=$7
outfile=$8
nucs_upstream=$9
ncpu=${10}

pwm_length_file=$outdir/pwm_length.txt

#----------------------------------------
# Data prep for final table
#----------------------------------------
# Compute PWM length
software/get_pwm_length.sh $pwms_file $pwm_length_file

# Reformat gene aliases
geneinfo_file=$outdir/genes_aliases.txt
awk -F '\t' '/^AT[0-9]G/{ if ($1 == gene) { s=s","$2; d=d","$3 } else { print gene"\t"s"\t"d; gene = $1; s=$2; d=$3} }' $genes_aliases_file | tail -n +2 | sed '1ilocus_name\tsymbol\tfull_name'  > $geneinfo_file

#----------------------------------------
# Parallel processing prep
#----------------------------------------
# split in features into smaller subsets

foutdir=$outdir/features
if [ ! -d $foutdir ]; then
	mkdir -p $foutdir
fi

header=header.txt
head -n 1 $features_rdat > $header
#header_line=`head -n 1 $features_rdat` 
#echo $header_line

nlines=100

prefix=$foutdir/feature.
#infile_noheader=t
#`tail -n +2 $features_rdat`

cat $features_rdat | tail -n +2 | split -d -l $nlines - $prefix

echo "Start to process tables $date"
count=1
for fin in `ls $prefix*`; do
	if [ $count -lt $ncpu ]; then
		# attach header line to each file
		#sed -i -e "1i$header_line"  $fin
		cat $header $fin > $fin.txt
		# convert the tables into long format - Each row has three items:(tss_id, factor_window, loglikScore)
		echo software/get_sites_table.R $fin.txt $geneinfo_file $peaks_file $pwm_length_file $tss_probs $model_weights 
		software/get_sites_table.R $fin.txt $geneinfo_file $peaks_file $pwm_length_file $tss_probs $model_weights $nucs_upstream $fin.long.txt &
                let "count=count+1"
        else
                echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
	fi
done
wait

echo "All tables are Done Processing!! $date"

tmpfile=$outdir/tmp.long.txt
# gather and unify all results
cat $prefix*.long.txt | grep -v "tss_id" > $tmpfile
head -n 1 $prefix"00.long.txt" | cat - $tmpfile > $outfile

echo "Cleaning up!"
#rm -f $prefix*
#rm -f -r $foutdir
rm -f $tmpfile





