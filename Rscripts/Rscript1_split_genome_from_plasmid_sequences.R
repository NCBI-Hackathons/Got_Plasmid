rm(list=ls())
getwd()



### Intputs ###
# ===========
my_gapless_files <- "gapless_genomes/fasta/gapless_genome_plasmid.fasta" #load file
#===========


### Code ###
#===========
row = scan(my_gapless_files,what=character(0),sep="\n")
chars = substr(row,1,1)
base = which(chars==">")
base = append(base, c(length(chars)))
  
  # Create a file with plasmid sequences
  plas <- base[grep("plasmid", row[base])] 

  for(i in 1:length(plas)){
   y <- which(plas[i] == base)
    genome <- row[base[y]:(base[y+1]-1)]
    gen <- chars[base[y]:(base[y+1]-1)]
    seq = gen!=">"
    seq = paste(genome[seq],collapse="")
    name <- gsub("Staphylococcus aureus subsp. aureus", "S_aureus",row[base[y]])
    name <- gsub(", .*", "",name)
    output = "plasmids/fasta/extra_plasmids.fasta"
    cat(name,sep="\n",file=output,append=TRUE)
    cat(seq,sep="\n",file=output,append=TRUE)
  }
  
  # megre into existing file: cat extra_plasmids.fasta >> all_plasmids.fasta

  # Create a file with gapless genome sequences
  gapless_gen <- base[grep("genome", row[base])]
  gapless_genomes <- gapless_gen[-grep("plasmid", row[gapless_gen])]
  for(i in 1:length(gapless_genomes)){
    y <- which(gapless_genomes[i] == base)
    genome <- row[base[y]:(base[y+1]-1)]
    gen <- chars[base[y]:(base[y+1]-1)]
    seq = gen!=">"
    seq = paste(genome[seq],collapse="")
    name <- gsub("Staphylococcus aureus subsp. aureus", "S_aureus",row[base[y]])
    name <- gsub(", .*", "",name)
    output = "gapless_genomes/fasta/gapless_genomes.fasta"
    cat(name,sep="\n",file=output,append=TRUE)
    cat(seq,sep="\n",file=output,append=TRUE)
  }
  
  # OPTIONAL: generate individual plasmid files from extra_plasmids.fasta
  

  extra.plas <- "plasmids/fasta/extra_plasmids.fasta"
  row = scan(extra.plas,what=character(0),sep="\n")
  chars = substr(row,1,1)
  base = which(chars==">")
  base = append(base, c(length(chars)+1))
  extra.plas.names <- as.data.frame(matrix(NA, (length(base)-1),1 ))

  
  for(i in 1:(length(base)-1)){
    seq <- row[base[i]:(base[i+1]-1)]
    name_plas <-  row[base[i]]
    name_plas <- gsub(">", "", name_plas)
    name_plas_refSeq <- gsub(" .*", "-", name_plas)
    name_plas <- paste("extra",gsub(" .*", "", name_plas), sep = "")
    extra.plas.names[i,1] <- name_plas_refSeq
    output = paste("plasmids/fasta/",name_plas,".fasta", sep = "")
    cat(seq,sep="\n",file=output,append=TRUE)
  }


  
#===========

  