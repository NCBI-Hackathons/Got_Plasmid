setwd("gapless_genomes/fasta/")

# makeblastdb_gapless_genomes.sh:
# -------------------------------


Script = "makeblastdb_gapless_genomes.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
cat(paste("makeblastdb -dbtype nucl -in gapless_genomes.fasta", sep = ""),sep="\n",file=Script,append=TRUE)
## chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid_gitHub/gapless_genomes/fasta/makeblastdb_gapless_genomes.sh
setwd("../../")