



Plasmid_ID <- scan("plasmids/assembly_ID/accession_plasmids.txt",what=character(0),sep="\t")
Plasmid_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "", s)})
Plasmid_fasta <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", ".fasta", s)})


# magicBlast one SRRA versus individual plasmid databases

setwd("plasmids/magic_output/")
Plasmid_ID <- scan("../assembly_ID/accession_plasmids.txt",what=character(0),sep="\t")
magicblast_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "_magic.txt", s)})

Script = "magicblast_plasmid.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for(i in 1:length(Plasmid_ID)){
  cat(paste("magicblast -db ../fasta/", Plasmid_fasta[i]," -sra SRR6227128 -no_unaligned -splice F -score 50 -outfmt tabular > ",magicblast_name[i] , sep =""),sep="\n",file=Script,append=TRUE)
}
## chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid_gitHub/plasmids/magic_output/magicblast_plasmid.sh
setwd("../../")
