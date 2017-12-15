

### Packages ###
#===============
library(data.table)
require(GenomicRanges)
library(IRanges)
#===============


### Functions ###
#===============
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}
#===============


### Inputs ###
#============


#--------- Extract gff information
setwd("../../gapless_genomes/gff/")
my_gff <- list.files(pattern = "\\.gff")
my_genome_names <- gsub("_genomic.gff","",my_gff)
my_genome_names <- gsub("1_.*","1",my_genome_names)
my_genome_names <- gsub("2_.*","2",my_genome_names)
my_genome_names <- gsub("3_.*","3",my_genome_names)

for (k in 1:length(my_gff)) {
  
  Gen.name <- my_gff[k]
  save.name <- paste(my_genome_names[k], "_Gff.csv", sep = "")
  
  #add columns of interest
  gff <- gffRead(gffFile = Gen.name)
  gff$GenBank_locus_tag <- getAttributeField(gff$attributes, "locus_tag")
  gff$Product <- getAttributeField(gff$attributes, "product")
  gff$Gene_name <- getAttributeField(gff$attributes, "gene")
  gff$Note <- getAttributeField(gff$attributes, "Note")
  gff$Pseudogene <- getAttributeField(gff$attributes, "gene_biotype")
  
  #remove rows and columns with non-relevant info
  gffClean <- gff[ ! gff$feature %in% c("CDS", "sequence_feature", "STS", "exon","transcript"), ]
  gffClean <- gffClean[-c(9,8,6,2)]
  
  gffClean$Pseudogene[which(gffClean$Pseudogene== "protein_coding")] <- NA
  gffClean$Pseudogene[which(gffClean$Pseudogene== "misc_RNA")] <- NA
  gffClean$Pseudogene[which(gffClean$Pseudogene== "rRNA")] <- NA
  gffClean$Pseudogene[which(gffClean$Pseudogene== "tRNA")] <- NA
  
  #add product info to main dataframe
  gffProduct <- gff[ gff$feature %in% c("CDS"), ]
  gffProduct <- gffProduct[c(which(colnames(gffProduct) == "start"), (which(colnames(gffProduct) == "Product")))]
  
  for (i in 1:nrow(gffProduct)){  
    if(length(gffProduct$start[i] >0)) {if (gffProduct$start[i] %in% gffClean$start ) {
      d <- which(gffProduct$start[i] == gffClean$start)
      gffClean$Product[d] <- gffProduct$Product[i]}
    }
  }
  
  #add gene size
  size.gene <- as.data.frame(matrix(data = NA, nrow = nrow(gffClean), ncol = 1))
  
  for (i in 1:nrow(gffClean)) {
    size.gene[1] <- gffClean$end - gffClean$start +1
  }
  
  colnames(size.gene) <- "Locus_size"
  gffClean <- cbind(gffClean, size.gene)
  x <- which(is.na(gffClean$Product) == T)
  gffClean$Product[x] <- gffClean$feature[x]
  genome_name <- gffClean$seqname[1]
  
  #save dataframe
  write.csv(gffClean, file = paste(genome_name,"-",save.name, sep = ""), row.names = FALSE, quote = FALSE) 
}

### Code ###
#============
setwd("../../")

SRA_names <- scan("SRA/SRA_ID.txt", what=character(0),sep="\n")
for (sr in 1:length(SRA_names)) {

setwd("plasmids/blast_vs_genome_output/")

SRA <- SRA_names[sr]

my_blast_Genome_files <- list.files(pattern = SRA)
my_blast_Genome_names <- gsub("_contigs_genomeBLAST.txt","",my_blast_Genome_files)
my_blast_Genome_list <- lapply(my_blast_Genome_files, function(i){fread( i, sep = "\t", header=T, data.table = F, fill = T, skip = 4)})

names(my_blast_Genome_list) <- my_blast_Genome_names
No_Genome_match <- my_blast_Genome_list[lapply(my_blast_Genome_list, nrow)==0]
Genome_match <- my_blast_Genome_list[lapply(my_blast_Genome_list, nrow)>0]

got_plas_table <- read.csv(paste("../../outputs/", SRA,"_got_plamid.csv", sep = ""))
#============

### Code ###
#===========
#--------- identify contigs with no match in gapless genomes
got_plas_table$genome_match <- NA
x <- which(got_plas_table$plasmid_ID %in% names(No_Genome_match))
got_plas_table$genome_match[x] <- "no_match"

#--------- identify contigs with matches in gapless genomes
column.names <- c("qseqid", "sseqid", "pident", "length_match", 
                  "mismatch", "gapopen","qstart", "qend", 
                  "sstart", "send", "evalue", "bitscore", "sgi", "stitle" )
Genome_match <- lapply(Genome_match, setNames, column.names)

#--------- Clean-up genome_match
for (i in 1: length(Genome_match)) {
      Genome_match[[i]] <- Genome_match[[i]][-grep("#", Genome_match[[i]]$qseqid),]
}


#--------- Add strand column and invert coordinates for "minus" strand
for (o in 1: length(Genome_match)) {
      Genome_match[[o]]$strand <- NA
      
      for (i in 1:nrow(Genome_match[[o]])) {
              if (Genome_match[[o]]$sstart[i] < Genome_match[[o]]$send[i]) {
              Genome_match[[o]]$strand[i] <- "plus"}else{ Genome_match[[o]]$strand[i] <- "minus"}
      }
          
      my_temp <- Genome_match[[o]]
      
      for (j in 1:nrow(Genome_match[[o]])) {
            if(my_temp$strand[j] == "minus"){
              my_temp$sstart[j] <- Genome_match[[o]]$send[j]
              my_temp$send[j] <- Genome_match[[o]]$sstart[j]}
      }
      Genome_match[[o]] <- my_temp
}

#--------- select best Bit score for each genome

Genome_match_temp1 <- list()
for (o in 1: length(Genome_match)) {
  uniq_contig <- unique(Genome_match[[o]]$qseqid)
  Genome_match_temp2 <- list()
      for (k in 1:length(uniq_contig)){
        d <- Genome_match[[o]][which(Genome_match[[o]]$qseqid == uniq_contig[k]),]
        uniq_gaplessGenome <- unique(d$sseqid)
        Genome_match_temp3 <- list()
            for (p in 1: length(uniq_gaplessGenome)){
              v <- d[which(d$sseqid == uniq_gaplessGenome[p]),]
              w <- v[order(v$bitscore, decreasing = T),]
              Genome_match_temp3[[p]] <- w[which(w$bitscore == max(w$bitscore)),]
            }
        Genome_match_temp3 <- do.call(rbind, Genome_match_temp3)
        Genome_match_temp2[[k]] <- Genome_match_temp3
      }
  Genome_match_temp2 <- do.call(rbind, Genome_match_temp2)
  Genome_match_temp1[[o]] <- Genome_match_temp2
}


names(Genome_match_temp1) <- names(Genome_match)
Genome_match <- Genome_match_temp1



#--------- identify affetced genes in genomes
for(o in 1:length(Genome_match)){
  setwd("../../gapless_genomes/gff/")
      Genome_match[[o]]$locus <- NA
      Genome_match[[o]]$gene <- NA
      Genome_match[[o]]$product <- NA
      Genome_match[[o]]$ref_gene_size <- NA
      Genome_match[[o]]$plasmid <- names(Genome_match)[o]
      
      uniq_genome <- unique(Genome_match[[o]]$stitle)
      uniq_genome_short <- sapply(uniq_genome, function(x) gsub("\\..*","",x))
  
      my_gff_files <- list.files(pattern = "\\_Gff.csv")
      my_gff <- sapply(my_gff_files, function(x) gsub("\\..*","",x))
  
  
      for(i in which(uniq_genome_short %in% my_gff)){
            nub <- which(Genome_match[[o]]$stitle == uniq_genome[i])
            my_gen <- Genome_match[[o]][nub,]
            f <- grep(uniq_genome_short[i], my_gff)
            
            my_gff_table <- read.csv(my_gff_files[f])
            my_gff_table$GenBank_locus_tag <- as.character(my_gff_table$GenBank_locus_tag)
            my_gff_table$Gene_name <- as.character(my_gff_table$Gene_name)
            my_gff_table$Product <- as.character(my_gff_table$Product)
            my_gff_table$Locus_size <- as.numeric(as.character(my_gff_table$Locus_size))

            
            ir <- IRanges(start = my_gen$sstart, end = my_gen$send)
            coord <- as.data.frame(reduce(ir))
            gff.ir <- IRanges(start = my_gff_table$start, end = my_gff_table$end)
            ol.gene <- findOverlaps(ir, gff.ir)
            ol.gene <- as.data.frame(ol.gene)
      
            for(j in 1:length(nub)){
              vec <- c("NA")
              vec <- append(vec,my_gff_table$GenBank_locus_tag[c(ol.gene$subjectHits[which(ol.gene$queryHits == j)])])
              vec <- vec[-1]
              vec <- paste(vec[!is.na(vec)], collapse = "/")
  
              my_gen$locus[j] <- vec
              
              vec2 <- c("NA")
              vec2 <- append(vec2,my_gff_table$Gene_name[c(ol.gene$subjectHits[which(ol.gene$queryHits == j)])])
              vec2 <- vec2[-1]
              vec2 <- paste(vec2[!is.na(vec2)], collapse = "/")

              my_gen$gene[j] <- vec2
              
              vec3 <- c("NA")
              vec3 <- append(vec3,my_gff_table$Product[c(ol.gene$subjectHits[which(ol.gene$queryHits == j)])])
              vec3 <- vec3[-1]
              vec3 <- paste(vec3[!is.na(vec3)], collapse = "/")

              my_gen$product[j] <- vec3
              
              vec4 <- c("NA")
              vec4 <- append(vec4,my_gff_table$Locus_size[c(ol.gene$subjectHits[which(ol.gene$queryHits == j)])])
              vec4 <- vec4[-1]
              vec4 <- paste(vec4[!is.na(vec4)], collapse = "/")

              my_gen$ref_gene_size[j] <- vec4
              
            }
            
            Genome_match[[o]][nub,] <- my_gen
      }
}



#--------- Clean-up to highlight intergenic regions
for(o in 1:length(Genome_match)){

      x <- which(Genome_match[[o]]$product == "region")
      Genome_match[[o]]$locus[x] <- "intergenic"
      Genome_match[[o]]$product[x] <- "-"
      Genome_match[[o]]$ref_gene_size[x] <- "-"
      Genome_match[[o]]$gene[x] <- "-"
      
      y <- which(is.na(Genome_match[[o]]$gene) == T)
      Genome_match[[o]]$gene[y] <- "not_available"
      y <- which(Genome_match[[o]]$gene == "")
      Genome_match[[o]]$gene[y] <- "not_available"
      
      y <- which(is.na(Genome_match[[o]]$locus) == T)
      Genome_match[[o]]$locus[y] <- "not_available"
      y <- which(Genome_match[[o]]$locus == "")
      Genome_match[[o]]$locus[y] <- "not_available"
      
      z <- which(Genome_match[[o]]$locus == Genome_match[[o]]$gene)
      Genome_match[[o]]$product[z] <- "not_available"
      Genome_match[[o]]$ref_gene_size[z] <- "not_available"
}

matching_genome_table <- do.call(rbind, Genome_match)
matching_genome_table <- matching_genome_table[c(ncol(matching_genome_table), 1:(ncol(matching_genome_table)-1))]

#--------- Save the table
write.csv(matching_genome_table,paste("../../outputs/",SRA,"_contig_vs_genomes.csv", sep = ""), row.names = F)
setwd("../../")
}
#============

 

