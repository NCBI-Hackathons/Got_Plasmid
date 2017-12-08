rm(list=ls())
getwd()
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/")

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
setwd("gff3_plasmid_files/")
my_gff <- list.files(pattern = "\\.gff")
my_genome_names <- gsub(".gff3","",my_gff)
#============


### Code ###
#===========

#--------- Extract gff information
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
    
        if(length(gffProduct$start[i] >0)) {
        
            if (gffProduct$start[i] %in% gffClean$start ) {
            d <- which(gffProduct$start[i] == gffClean$start)
            gffClean$Product[d] <- gffProduct$Product[i]
            }
        
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

#============



