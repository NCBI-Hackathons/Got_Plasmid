rm(list=ls())


#### from https://github.com/gschofl/biofiles

setwd("plasmids/")

### Packages ###
#===============

library(devtools)
library(biofiles)
library(reutils)
library(IRanges)
require(GenomicRanges)
#===============


### Intputs ###
# =============


gbk.acc <- "assembly_ID/accession_plasmids.txt"
#dowloaded from NCBI Nucleotide using:
# plasmid[title] AND staphylococcus[title]
# Filters: Species = bacteria
# Filters: Molecule types = genomic RNA/DNA
# Filters: Genetic compartments = Plasmid
# download accession table
#=============


### Code ###
# =============
acc_files <- read.delim(gbk.acc, sep = "\t", header = F, stringsAsFactors = F)
#acc_files <- names(output_list)[1:10]
save.list <- list()

#for (i in 1:nrow(acc_files)) {
for (i in 1:nrow(acc_files)) {
  
#acc <- as.character(acc_files[i,1])
acc <- acc_files[i,1]

gb_file <- reutils::efetch(acc, "nuccore", rettype = "gbwithparts", retmode = "text")
rec <- biofiles::gbRecord(gb_file)
f <- biofiles::ft(rec)

source <- biofiles::filter(rec, key = "source")
plasmid_ranges <- biofiles::ranges(source,include = c("locus_tag", "gene", "protein_id", "product", "note","old_locus_tag","pseudo"))

genes <- biofiles::filter(rec, key = "gene")
genes_ranges <- biofiles::ranges(genes, include = c("locus_tag", "gene", "protein_id", "product", "note","old_locus_tag","pseudo"))

CDS <- biofiles::filter(rec, key = "CDS")
CDS_ranges <- biofiles::ranges(CDS, include = c("locus_tag", "gene", "protein_id", "product", "note","old_locus_tag","pseudo"))

tRNA <- biofiles::filter(rec, key = "tRNA")
tRNA_ranges <- biofiles::ranges(tRNA, include = c("locus_tag", "gene", "protein_id", "product", "note","old_locus_tag","pseudo"))

df0 <- as.data.frame(plasmid_ranges, row.names = NULL)
df1 <- as.data.frame(genes_ranges, row.names = NULL)
df2 <- as.data.frame(CDS_ranges,row.names = NULL)
df3 <- as.data.frame(tRNA_ranges,row.names = NULL)

if(length(which(duplicated(names(genes_ranges[1:length(genes_ranges)])))) != 0){
  x <- which(duplicated(names(genes_ranges[1:length(genes_ranges)])))
  for(fs in 1:length(x)){
    names(genes_ranges)[x[fs]] <- paste(names(genes_ranges)[x[fs]], "_duplicate",sep="")
  }
}

if(length(which(duplicated(names(CDS_ranges[1:length(CDS_ranges)])))) != 0){
  x <- which(duplicated(names(CDS_ranges[1:length(CDS_ranges)])))
  for(fs in 1:length(x)){
    names(CDS_ranges)[x[fs]] <- paste(names(CDS_ranges)[x[fs]], "_duplicate",sep="")
  }
}
      
if(ncol(df0) != ncol(df1)){
  df1 <- as.data.frame(genes_ranges)
  df1$locus_tag <- rownames(df1)
  df1$gene <- NA
  df1$protein_id <- NA
  df1$product <- NA
  df1$note <- NA
  df1$old_locus_tag <- NA
  df1$pseudo <- NA
}

if(ncol(df0) != ncol(df2)){
  df2 <- as.data.frame(CDS_ranges)
  df2$locus_tag <- rownames(df1)
  df2$gene <- NA
  df2$protein_id <- NA
  df2$product <- NA
  df2$note <- NA
  df2$old_locus_tag <- NA
  df2$pseudo <- NA
}

  
  

df <- rbind(df0, df1,df2)
if(nrow(df3)>0){
  df <- rbind(df, df3)
}
  save.list[[i]] <- df

}

GenBank.table <- do.call(rbind,save.list)
GenBank.table$start <- as.numeric(GenBank.table$start)
GenBank.table$end <- as.numeric(GenBank.table$end)
GenBank.table$width <- as.numeric(GenBank.table$width)
GenBank.table$temp <- paste(GenBank.table$seqnames, GenBank.table$end)


#add product info to main dataframe
gbkProduct <- GenBank.table[ GenBank.table$key %in% c("CDS"), ]

for (i in 1:nrow(gbkProduct)){  

  if (gbkProduct$temp[i] %in% GenBank.table$temp ) {
    d <- which(gbkProduct$temp[i] == GenBank.table$temp)
    GenBank.table$product[d] <- gbkProduct$product[i]
    GenBank.table$note[d] <- gbkProduct$note[i]
    GenBank.table$old_locus_tag[d] <- gbkProduct$old_locus_tag[i]
    GenBank.table$pseudo[d] <- gbkProduct$pseudo[i]
    GenBank.table$protein_id[d] <- gbkProduct$protein_id[i]
  }
  
}

gbkClean <- GenBank.table[ ! GenBank.table$key %in% c("CDS"), ]

#add missing CDS

uniq.genome <- unique(GenBank.table$seqnames)

for (i in 1:length(uniq.genome)) {
  w <- gbkProduct$end[which(gbkProduct$seqnames == uniq.genome[i])]
  v <- gbkClean$end[which(gbkClean$seqnames == uniq.genome[i])]
  z <- w %in% v
  x <- which(gbkProduct$seqnames == uniq.genome[i])
  
  gbkClean <- rbind(gbkClean, gbkProduct[x[which(z == F)],])
  
}

write.table(gbkClean, file = "genBank/all_plasmids_GenBank.txt", row.names = FALSE, quote = FALSE, sep = "\t")
setwd("../")
# =============