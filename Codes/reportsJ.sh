#!/bin/bash

#######################################################
# Parameters meaning
# gspmax: max. number of gapped HSPs (GSPs) saved per subject sequence (default 0; 0 => unlimited)
# -B -V: the B and V options limit the number of subject sequences for which any results whatsoever are reported, regardless of the number of HSPs or GSPs found. 
# -E: Expectation value (E) threshold for saving hits 
# -cpus: no. of processors to utilize on multi-processor systems
#######################################################


# Search protein data (Uniprot database)
# If you decided to download your own datasets, you need to make database for the blastx
xdformat -p -k uniprot_sprot.fasta
blastx ./report_run/sprot notKnown.fa  -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.spwb
python ./report_run/wublastx2gff.py notKnown.fa.spwb > notKnown.fa.spwb.gff

# Search from GB_TE database
# If you decided to download your own datasets, you need to make database for the blastx
xdformat -p -k GB_TE.21032016.fa -o GB_TE.new
blastx ./BlastDB/GB_TE.new notKnown.fa -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.tewb
python ./report_run/wublastx2gff.py notKnown.fa.tewb > notKnown.fa.tewb.gff

# Search Retrovirus data
# If you decided to download your own datasets, you need to make database for the tblastx
xdformat -p -k all_retrovirus.fasta
tblastx ./BlastDB/all_retrovirus.fasta notKnown.fa -gspmax=1 -E 0.00001 -B 1 -V 1 -cpus=32 > notKnown.fa.ervwb
python ./report_run/wublastx2gff.py notKnown.fa.ervwb > notKnown.fa.ervwb.gff
