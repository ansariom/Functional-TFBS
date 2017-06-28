#!/usr/bin/perl
#This script reads wiggle file and extracts the probability for TSS location (one location)

$wigFile = shift;
$peakName = shift;

open(IN, $wigFile) || die "Can't open $wigFile";

while (<IN>) {	
	$line = $_;
	chomp $line;
	if ($line =~ m/^0/g) {
	     $prob = $line;
	      if ($prob > 0 ){
	          print $peakName . "\t". $prob . "\n";
	      }
        }
}
close(IN);
