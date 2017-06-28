#!/bin/bash

pwm_file=$1
outfile=$2

awk '{if ($0 ~ /^>/) {pwm = substr($0,3); count = 0} else if ($0 ~ /^=/){ count=count+1} else {print pwm"\t"count; count=0} }' $pwm_file > $outfile
