#!/bin/bash

for i in `cat ../assembly_ID/S_aureus_gapless_assemblies`; do wget `esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | 
awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`; done 

gunzip *.gz

cat *.fna > S_aureus_gapless_genome_plasmid.fasta

mv S_aureus_gapless_genome_plasmid.fasta ../genome_fasta/.

rm *.fna

