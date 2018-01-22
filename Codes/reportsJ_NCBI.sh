#!/bin/bash

#######################################################
# Parameters meaning
# -db: BLAST database name
# -query: Input file name
# -max_hsps: Set maximum number of HSPs per subject sequence to save for each query, ncbi-blast doesn't have gsps option, but we tested with hsps in wu-blast, the result remain almost the same
# -seg Filter query sequence with SEG (Format: 'yes', 'window locut hicut', or 'no' to disable), default wu-blast is off
# -evalue: Expectation value (E) threshold for saving hits 
# -num_threads: Number of threads (CPUs) to use in the BLAST search
# -max_target_seqs: Maximum number of aligned sequences to keep 
# -word_size: Word size for wordfinder algorithm. The wu-blast default is 3, we've tried both word_size 2 and 3 in NCBI-blast, 2 can find more same results compared to wublast
########################################################


# Search protein data (Uniprot database)
# If you decided to download your own datasets, you need to make database for the blastx
makeblastdb -in uniprot_sprot.fasta -dbtype nucl
blastx -db uniprot_sprot.fasta -query notKnown.fa -max_hsps 1 -seg no -evalue 0.00001 -num_threads 32 -max_target_seqs 1 -word_size 2 -outfmt 6 -out notKnown.fa.spwb.ncbi
awk '{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target sp|"$2" "$9" "$10}' notKnown.fa.spwb.ncbi > tmp
awk '{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}' tmp > notKnown.fa.spwb.gff

# Search from GB_TE databae
# If you decided to download your own datasets, you need to make database for the blastx
makeblastdb -in GB_TE.21032016.fa -dbtype prot -out GB_TE.new
tblastx -db GB_TE.new -query notKnown.fa -max_hsps 1 -seg no -evalue 0.001 -num_threads 32 -max_target_seqs 1 -word_size 2 -outfmt 6 -out notKnown.fa.tewb.ncbi
awk '{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target sp|"$2" "$9" "$10}' notKnown.fa.ervwb.ncbi > tmp
awk '{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}' tmp > notKnown.fa.tewb.gff
# Search Retrovirus data
# If you decided to download your own datasets, you need to make database for the tblastx
makeblastdb -in all_retrovirus.fasta -dbtype nucl 
tblastx -db all_retrovirus.fasta -query notKnown.fa -max_hsps 1 -seg no -evalue 0.001 -num_threads 32 -max_target_seqs 1 -word_size 2 -outfmt 6 -out notKnown.fa.ervwb.ncbi
awk '{print $1"\t""blast""\t""hit""\t"$7"\t"$8"\t"$11"\t"".""\t"".""\t""Target sp|"$2" "$9" "$10}' notKnown.fa.ervwb.ncbi > tmp
awk '{if($4>$5) print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9" "$10" "$11" "$12; else print $0}' tmp > notKnown.fa.ervwb.gff