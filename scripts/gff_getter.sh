#!/bin/bash
for i in cat 'S_aureus_gapless_assemblies'; do wget 'esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.gff.gz"}'; done
gunzip *.gz
