#!/bin/bash
for i in `cat gapless_genome_assemblies.txt`; do wget `esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`; done
gunzip *.gz
cat *.fna > gapless_genome_plasmid.fasta
mv gapless_genome_plasmid.fasta ../fasta/.
rm *.fna
