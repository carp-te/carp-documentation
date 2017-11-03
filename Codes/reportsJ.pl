#!/usr/bin/perl -w

use strict;
my $dataDir = "./";
my $ECrit = 0.00001;

#print "Enter name of file listing libraries to process: ";
my $libs = '/home/a1635743/libs.txt';

open (INFILE, "$libs") || die "You dope, no input file by that name, try again...";
 
while (<INFILE>) {
  chomp;
  system ("blastx /home/a1635743/report_run/sprot $dataDir$_  -gspmax=1 -E $ECrit -B 1 -V 1 -cpus=32 > $_.spwb");
  system ("python /home/a1635743/report_run/wublastx2gff.py $_.spwb > $_.spwb.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.spwb.gff > $_.spwb.txt");
  system ("/home/a1635743/report_run/piler -annot $_.gff -rep $_.spwb.gff -out $_.spwbannot.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.spwbannot.gff > $_.spwbannot.txt");

  system ("blastx /home/a1635743/BlastDB/GB_TE.new $dataDir$_ -gspmax=1 -E $ECrit -B 1 -V 1 -cpus=32 > $_.tewb");
  system ("python /home/a1635743/report_run/wublastx2gff.py $_.tewb > $_.tewb.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.tewb.gff > $_.tewb.txt");
  system ("/home/a1635743/report_run/piler -annot $dataDir$_.gff -rep $_.tewb.gff -out $_.tewbannot.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.tewbannot.gff > $_.tewbannot.txt");

  system ("tblastx /home/a1635743/BlastDB/all_retrovirus.fasta $dataDir$_ -gspmax=1 -E $ECrit -B 1 -V 1 -cpus=32 > $_.ervwb");
  system ("python /home/a1635743/report_run/wublastx2gff.py $_.ervwb > $_.ervwb.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.ervwb.gff > $_.ervwb.txt");
  system ("/home/a1635743/report_run/piler -annot $dataDir$_.gff -rep $_.ervwb.gff -out $_.ervwbannot.gff");
  system ("python /home/a1635743/report_run/gff_pretty.py $_.ervwbannot.gff > $_.ervwbannot.txt");

  chdir ("..");
	}
