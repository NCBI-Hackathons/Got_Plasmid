setwd("plasmids/contigs/")
contigs_fasta <- list.files(pattern = "\\_matching_contigs.fasta")
contigs_fasta <- contigs_fasta[-grep("fasta.", contigs_fasta)]
contigs_blast <- gsub("_matching_contigs.fasta","_contigs_genomeBLAST.txt",contigs_fasta)

Script = "blast_contigs_vs_gapless_genomes.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for (i in 1: length(contigs_fasta)) {
  cat(paste("blastn -db ../../gapless_genomes/fasta/gapless_genomes.fasta -num_threads 1 -max_target_seqs 5 -outfmt '7 std sgi stitle' -evalue 0.0000000001 -query ", contigs_fasta[i]," -out ../blast_vs_genome_output/", contigs_blast[i], sep = ""),sep="\n",file=Script,append=TRUE)
}

## chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid/scripts/blast_contigs_vs_gapless_genomes.sh
setwd("../../")