setwd("plasmids/assembly_ID/")
Plasmid_ID <- scan("accession_plasmids.txt",what=character(0),sep="\t")
Plasmid_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "", s)})
Plasmid_fasta <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", ".fasta", s)})

# esearch_plasmid.sh:
# -------------------
Script = "esearch_plasmid.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(Plasmid_ID)){
  cat(paste("esearch -db nucleotide -query ", Plasmid_ID[i]," | efetch -format fasta > ",Plasmid_fasta[i], sep = ""),sep="\n",file=Script,append=TRUE)
}
cat(paste("cat *.fasta > S_aureus_all_plasmids.fasta"), sep="\n",file=Script,append=TRUE)
cat(paste("mv *.fasta ../fasta/."), sep="\n",file=Script,append=TRUE)
## chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid_gitHub/plasmids/assembly_ID/esearch_plasmid.sh
setwd("../../")
