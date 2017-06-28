#!/bin/bash

# Modified by Mitra - this is the most up to date version
#---------------------------------------------------------
# 1- Old Java cals is been replaced by tfbs_scan.jar
# 2- Input parameters are changed and named accordingly


########################################################################################
# * Copyright (c) 2014 Oregon State University
# * Authors: Taj Morton
# *
# * This program is distributed under the terms listed in the
# * LICENSE file included with this software. 
# *
# * IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT,
# * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
# * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF OREGON
# * STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. OREGON STATE
# * UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# * AND ANY STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER
# * IS ON AN "AS IS" BASIS, AND OREGON STATE UNIVERSITY HAS NO OBLIGATIONS TO
# * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
# *
# * Contact: megrawm@science.oregonstate.edu
########################################################################################
modelFilename=$1  # Model file is a sparse matrix containing model coefficients
seqsFile=$2     # Sequence to scan 
locsFile=$3     # A bed file containing locations of each sequence
pwm_file=$4	  # PWM File
roeFwd=$5	  # ROE fed
roeRev=$6	  # ROW rev
scan_hw=$7	# half width of scan region centered at TSS (0)
promoter_hw=$8	# needs to compute the genomic location that wig starts
outputDir=$9	  # Where to place results
scripts_dir=${10}
l1logreg_path=${11}


let "wig_start_loc = $promoter_hw - $scan_hw"
wig_step=1
bgWin=250
smoothHalfWin=2

seqsBasename="`basename $seqsFile`"
scansOut="$outputDir/$seqsBasename.features"
resultsDir="$outputDir/classification/$seqsBasename"
scriptOutputDir="$outputDir/classifier_stdout/$seqsBasename"
wigDir="$outputDir/scans/$seqsBasename"

function rdatToMM() {
    infile="$1"
    outfile="$2"
    $scripts_dir/Rdat_to_sparse_test_noclass.pl "$infile" "$outfile"
}

function generateScans() {
    outLocation=`readlink -f $1`
    for seqName in `grep '^>' $seqsFile |sed -e 's/^>\([^[:space:]]\+\).*$/\1/'`; do        
	#cd "$LOGLIK_PATH/classes"
        #java loglikscan.SumScoreVarsSeqLocalBG2FastWVEF_MP "Range 1 -$nucsScanFromTSS $nucsScanFromTSS" $downstream_len $roeFwd $roeRev $thisDir/$seqsFile $pwm_zeroes_file $pwm_file $bgWin $outLocation/$seqName.Rdat "$seqName"
        java -jar software/tfbs_scan.jar GenFeatures --pos "Range 1 -$scan_hw $scan_hw" -B $bgWin -N $promoter_hw --seqName "$seqName" $roeFwd $roeRev $seqsFile $pwm_file $outLocation/$seqName.Rdat
	

#        cd "$thisDir"
        rdatToMM "$outLocation/$seqName.Rdat" "$outLocation/$seqName.mm"
    done
}

function classifyExample() {
    modelFilename="$1"
    featuresFile="$2"
    resultsFile="$3"
    outputFile="$4"

    echo l1_logreg_classify -p "$modelFilename" "$featuresFile" "$resultsFile" 2>&1 | tee "$outputFile"
    $l1logreg_path/l1_logreg_classify -p "$modelFilename" "$featuresFile" "$resultsFile" 2>&1 | tee "$outputFile"
}

function smoothScan() {
    resultsFile="$1"
    outputFile="$2"
    smoothWin="$3"

    $scripts_dir/results_to_smoothed_results.pl "$smoothWin" "$resultsFile" "$outputFile"
}

function resultsToWig() {
    resultsFile="$1"
    locsFile="$2"
    wigDir="$3"
    $scripts_dir/results_to_wig3.pl $wig_start_loc $wig_step "$resultsFile" "$locsFile"  "$wigDir"
}

mkdir -p $scansOut
generateScans $scansOut

mkdir -p "$resultsDir"
mkdir -p "$outputDir"
mkdir -p "$scriptOutputDir"
mkdir -p "$wigDir"
for peak in $scansOut/*.mm; do
    peak_name=`basename $peak`
    classifyExample "$modelFilename" "$peak" "$resultsDir/$peak_name.result" "$scriptOutputDir/$peak_name.output"

    if [ ! -f "$resultsDir/$peak_name.result" ]; then
        echo "Classification failed."
        echo "Is l1_logreg_classify installed, and can it be called from the current directory?"
        echo "Check $outputDir/$peak_name.output for more error information."
    else
        smoothScan "$resultsDir/$peak_name.result" "$resultsDir/$peak_name.result.smoothed" "$smoothHalfWin" 
        resultsToWig "$resultsDir/$peak_name.result.smoothed" "$locsFile" "$wigDir"
    fi
done

