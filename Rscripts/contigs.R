setwd("plasmids/contigs/")
contig_files = list.files(pattern = "\\_matching_contigs.fasta")
contig_names = sapply(contig_files, 
                      function(s) {s = gsub("_matching_contigs.fasta", "", s)})

# make individual contig db:
# ----------------------------------
Formatblast = "Formatblast_contigs.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Formatblast,append=TRUE)
for(i in 1:length(contig_names)) {
  cat(paste("makeblastdb -dbtype nucl -in ", contig_files[i], sep=""),sep="\n",file=Formatblast,append=TRUE)
}

# Blast contigs against contig databases:
# --------------------------------------
Script = "reciprocal_contig_Blast.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(contig_names)) {
  for(j in 1:length(contig_names))  {
    cat(paste("blastn -db ", contig_files[i], " -num_threads 24 -max_target_seqs 5 -outfmt '6 std sgi stitle' -evalue 0.0000000001 -query  ",contig_files[j], " -out reciprocal_blast/",contig_names[j],"vs",contig_names[i],".txt", sep=""),sep="\n",file=Script,append=TRUE)
    
  }
}
# chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid_gitHub/plasmids/contigs/*.sh
setwd("../../")