
Plasmid_ID <- scan("plasmids/assembly_ID/accession_plasmids.txt",what=character(0),sep="\t")
Plasmid_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "", s)})
Plasmid_fasta <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", ".fasta", s)})

# esearch_plasmid.sh:


setwd("plasmids/fasta/")

# makeblastdb_individual_plasmid.sh:
# ----------------------------------
Script = "makeblastdb_individual_plasmid.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(Plasmid_ID)){
  cat(paste("makeblastdb -dbtype nucl -in ", Plasmid_fasta[i], sep = ""),sep="\n",file=Script,append=TRUE)
}
setwd("../../")