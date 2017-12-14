#!/bin/bash
for i in $(cat gapless_genome_assemblies.txt); do wget $(esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.gff.gz"}'); done
gunzip *.gz
mv *.gff ../gff/.
