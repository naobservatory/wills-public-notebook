#!/usr/bin/env bash

# Specify input and output paths
infile=$1
outfile=$2

# Run blastn
blastn -query $1 -out $2 -db nt -remote -outfmt "6 qseqid sseqid sgi staxid qlen evalue bitscore qcovs length pident mismatch gapopen sstrand qstart qend sstart send" -perc_identity 60 -max_hsps 5 -num_alignments 250 -qcov_hsp_perc 30
