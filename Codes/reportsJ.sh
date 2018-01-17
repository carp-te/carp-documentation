#!/bin/bash

# Filter protein data (Uniprot database)
blastx ./report_run/sprot notKnown.fa  -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.spwb
python ./report_run/wublastx2gff.py notKnown.fa.spwb > notKnown.fa.spwb.gff
python ./report_run/gff_pretty.py notKnown.fa.spwb.gff > notKnown.fa.spwb.txt
report_./run/piler -annot notKnown.fa.gff -rep notKnown.fa.spwb.gff -out notKnown.fa.spwbannot.gff
python ./report_run/gff_pretty.py notKnown.fa.spwbannot.gff > notKnown.fa.spwbannot.txt

# Identify from GB_TE databae
blastx ./BlastDB/GB_TE.new notKnown.fa -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.tewb
python ./report_run/wublastx2gff.py notKnown.fa.tewb > notKnown.fa.tewb.gff
python ./report_run/gff_pretty.py notKnown.fa.tewb.gff > notKnown.fa.tewb.txt
./report_run/piler -annot notKnown.fa.gff -rep notKnown.fa.tewb.gff -out notKnown.fa.tewbannot.gff
python ./report_run/gff_pretty.py notKnown.fa.tewbannot.gff > notKnown.fa.tewbannot.txt

# Filter Retrovirus data
tblastx ./BlastDB/all_retrovirus.fasta notKnown.fa -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.ervwb
python ./report_run/wublastx2gff.py notKnown.fa.ervwb > notKnown.fa.ervwb.gff
python ./report_run/gff_pretty.py notKnown.fa.ervwb.gff > notKnown.fa.ervwb.txt
./report_run/piler -annot notKnown.fa.gff -rep notKnown.fa.ervwb.gff -out notKnown.fa.ervwbannot.gff
python ./report_run/gff_pretty.py notKnown.fa.ervwbannot.gff > notKnown.fa.ervwbannot.txt