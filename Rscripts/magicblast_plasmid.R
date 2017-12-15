
SRA <- scan("SRA/SRA_ID.txt", what=character(0),sep="\n")
Plasmid_ID <- scan("plasmids/assembly_ID/accession_plasmids.txt",what=character(0),sep="\t")
Plasmid_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "", s)})
Plasmid_fasta <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", ".fasta", s)})

      # magicBlast one SRRA versus individual plasmid databases
      setwd("plasmids/magic_output/")
      magicblast_name <- sapply(Plasmid_ID, function(s) {s = gsub("\\..*", "_magic.txt", s)})
      
      Script = "magicblast_plasmid.sh"
      cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
      for(j in 1:length(SRA)){
          for(i in 1:length(Plasmid_ID)){
            cat(paste("magicblast -db ../fasta/", Plasmid_fasta[i]," -sra ", SRA[j]," -no_unaligned -splice F -score 50 -outfmt tabular > ",SRA[j], "-", magicblast_name[i] , sep =""),sep="\n",file=Script,append=TRUE)
          }
      }

setwd("../../")

