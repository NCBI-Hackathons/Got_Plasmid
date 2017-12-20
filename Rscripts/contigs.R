

SRA <- scan("SRA/SRA_ID.txt", what=character(0),sep="\n")
setwd("plasmids/contigs/")
contig_files = list.files(pattern = "\\_matching_contigs.fasta")
contigs_fasta <- contig_files[-grep("fasta.", contig_files)]
contig_names = sapply(contigs_fasta, 
                      function(s) {s = gsub("_matching_contigs.fasta", "", s)})

# make individual contig db:
# ----------------------------------
Formatblast = "Formatblast_contigs.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Formatblast,append=TRUE)
for(i in 1:length(contig_names)) {
  cat(paste("makeblastdb -dbtype nucl -in ", contigs_fasta[i], sep=""),sep="\n",file=Formatblast,append=TRUE)
}

# Blast contigs against contig databases:
# --------------------------------------
Script = "reciprocal_contig_Blast.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(contig_names)) {
  for(j in 1:length(contig_names))  {
    cat(paste("blastn -db ", contigs_fasta[i], " -num_threads 1 -max_target_seqs 5 -outfmt '6 std sgi stitle' -evalue 0.0000000001 -query  ",contigs_fasta[j], " -out reciprocal_blast/",contig_names[j],"vs",contig_names[i],".txt", sep=""),sep="\n",file=Script,append=TRUE)
    
  }
}
setwd("../../")