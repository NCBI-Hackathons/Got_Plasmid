setwd("plasmids/fasta/")

# makeblastdb_individual_plasmid.sh:
# ----------------------------------
Script = "makeblastdb_individual_plasmid.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(Plasmid_ID)){
  cat(paste("makeblastdb -dbtype nucl -in ", Plasmid_fasta[i], sep = ""),sep="\n",file=Script,append=TRUE)
}
setwd("../../")