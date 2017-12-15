
### Packages ###
#===============
library(data.table)
require(GenomicRanges)
library(IRanges)
#===============


### Functions ###
#===============
readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n")
  chars = substr(row,1,1)
  base = chars!=">"
  seq = paste(row[base],collapse="")
  return(toupper(unlist(strsplit(seq,""))))
}
#===============


### Inputs ###
#============
SRA_names <- scan("SRA/SRA_ID.txt", what=character(0),sep="\n")


for (sr in 1:length(SRA_names)) {
setwd("plasmids/magic_output/")
SRA = SRA_names[sr]

my_magic_files <- list.files(pattern = SRA)
my_magic_names <- gsub(".magic.txt","",my_magic_files)
plasmid_names <- gsub(".*-","",my_magic_names)
my_magic_list <- lapply(my_magic_files, function(i){fread( i, sep = "\t", header=T, data.table = F, fill = T, skip = 2)})

names(my_magic_list) <- my_magic_names
my_magic_list <- my_magic_list[lapply(my_magic_list, nrow)>0]
#============

### Code ###
#===========
output_list <- list()

#--------- loop over each magicBlast results
for (o in 1:length(my_magic_list)) {
  
  my_magic_table <- my_magic_list[[o]]
  my_magic_plas_name <- plasmid_names[o]
  ref <- readFastaRef(paste("../fasta/", my_magic_plas_name, ".fasta", sep = ""))
  lgth_ref <- length(ref)  
  
  #--------- Add strand columns
  my_magic_table$contigs_number <- NA
  my_magic_table$overlapping_contigs <- NA
  my_magic_table$overl_coord_start <- NA
  my_magic_table$overl_coord_end <- NA
  my_magic_table$overl_size <- NA
  my_magic_table$ref_size <- lgth_ref
  my_magic_table$overall_coverage <- NA
  my_magic_table <- my_magic_table[order(my_magic_table$`reference start`),]
  
  #--------- invert coordinates of "-" strand
  my_temp <- my_magic_table
  
  for (i in 1:nrow(my_magic_table)) {
    
    if(my_temp$`reference strand`[i] == "minus"){
      my_temp$`reference start`[i] <- my_magic_table$`reference end`[i]
      my_temp$`reference end`[i] <- my_magic_table$`reference start`[i]
    }
    
  }
  
  my_magic_table <- my_temp
  
  #--------- calculate coordinates of overlapping reads
  ir <- IRanges(start = my_magic_table$`reference start`, end = my_magic_table$`reference end`)
  ol <- findOverlaps(ir, reduce(ir))
  ol <- as.matrix(ol)
  coord <- as.data.frame(reduce(ir))
  
  my_magic_table$overlapping_contigs <- ol[,2]
  my_magic_table$contigs_number <- paste("contig_",ol[,2], sep = "")
  
  for (i in 1:nrow(coord)) {
    my_magic_table$overl_coord_start[which(my_magic_table$overlapping_contigs == i)] <- coord[i,1]
    my_magic_table$overl_coord_end[which(my_magic_table$overlapping_contigs == i)] <- coord[i,2]
    my_magic_table$overl_size[which(my_magic_table$overlapping_contigs == i)] <- coord[i,3]
  }
  
  my_magic_table$overall_coverage <- round(sum(as.numeric(coord[,3]))/lgth_ref*100)
  my_magic_table <- my_magic_table[-c(2,4,5,6,11,14,18,19,23,24)]
  output_list[[o]] <- my_magic_table
  names(output_list)[o] <- my_magic_plas_name
  
}

#============

setwd("../../")

### Cricos ###
#============
#--------- load gbk


my_gbk <- read.delim("plasmids/genBank/all_plasmids_GenBank.txt", header = T, sep = "\t")
my_gbk$start <- as.numeric(as.character(my_gbk$start))
my_gbk$end <- as.numeric(as.character(my_gbk$end))
my_gbk$width <- as.numeric(as.character(my_gbk$width))
my_gbk$simplify_genome_name <- sapply(my_gbk$seqnames, function(x) gsub(" .*","",x))



x <- which(names(output_list) %in% my_gbk$simplify_genome_name)


for (r in 1:length(x)) {
  
      nam <- names(output_list[x[r]])
      col.start <- which(colnames(output_list[[x[r]]]) == "contigs_number")
      col.end <- ncol(output_list[[x[r]]])
      new.df <- output_list[[x[r]]][, col.start:col.end]
      nodupl.df <- new.df[!duplicated(new.df[,1]),]
      my_gbk_table <- my_gbk[which(my_gbk$simplify_genome_name == nam),]
      size.plas <- max(my_gbk_table$end)
      
      #--------- create Karyotype
      setwd("circos/circos-0.69-6/circos_plasmid/data/")
      karyo <- paste(SRA,"-",nam, "_karyotpye.txt", sep = "")
      cat(paste("chr - chr1 1 0 ", size.plas, " black",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 0 ", round(size.plas/10*1), " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*1), " ", round(size.plas/10*2), " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*2), " ", round(size.plas/10*3), " lgrey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*3), " ", round(size.plas/10*4), " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*4), " ", round(size.plas/10*5), " lgrey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*5), " ", round(size.plas/10*6), " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*6), " ", round(size.plas/10*7), " lgrey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*7), " ", round(size.plas/10*8), " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*8), " ", round(size.plas/10*9), " lgrey",sep=""),sep="\n",file = karyo,append=TRUE)
      cat(paste("band chr1 1.1 1.1 ", round(size.plas/10*9), " ", size.plas, " grey",sep=""),sep="\n",file = karyo,append=TRUE)
      
      #--------- create contig highlights.txt
      highlight1 <- paste(SRA,"-",nam, "_highlight1.txt", sep = "")
            for (i in 1:round(nrow(nodupl.df))) {
                
                cat(paste("chr1", nodupl.df$overl_coord_start[i], 
                          nodupl.df$overl_coord_end[i],
                          paste("fill_color=chr", sample(1:20,1), sep=""), 
                          sep="\t"),
                    sep="\n",file = highlight1,append=T)
            }
      
      #--------- create cover.txt
      cover1 <- paste(SRA,"-",nam, "_cover1.txt", sep = "")
            for (i in 1:round(nrow(nodupl.df))) {
              
              cat(paste("chr1", nodupl.df$overl_coord_start[i], 
                        nodupl.df$overl_coord_end[i],
                        paste("fill_color=vlgrey",sep=""), 
                        sep="\t"),
                  sep="\n",file = cover1,append=T)
            }
      
      #--------- create gene_block.txt
      gene.blocks <- paste(SRA,"-",nam, "_gene_block.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        paste("fill_color=chr", sample(1:20,1), sep=""),
                        sep="\t"),
                  sep="\n",file = gene.blocks,append=T)
            }
      temp <- read.delim(paste(SRA,"-",nam, "_gene_block.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      write.table(temp, paste(SRA,"-",nam, "_gene_block.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create label_locus.txt
      label.locus <- paste(SRA,"-",nam, "_locus_labels.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", 
                        my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        my_gbk_table$locus_tag[i],
                        sep="\t"),
                  sep="\n",file = label.locus,append=T)
            }
      temp <- read.delim(paste(SRA,"-",nam, "_locus_labels.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      write.table(temp, paste(SRA,"-",nam, "_locus_labels.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create label_genes.txt
      label.genes <- paste(SRA,"-",nam, "_gene_labels.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", 
                        my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        my_gbk_table$gene[i],
                        sep="\t"),
                  sep="\n",file = label.genes,append=T)
            }
      temp <- read.delim(paste(SRA,"-",nam, "_gene_labels.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      s <- which(is.na(temp$V4) == T)
            if (length(s)>0) {
              temp <- temp[-which(is.na(temp$V4) == T),]
            }
      write.table(temp, paste(SRA,"-",nam, "_gene_labels.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create separator.txt
      separator <- paste(SRA,"-",nam, "_separator.txt", sep = "")
        cat(paste("chr1", 
                  0, 
                  max(my_gbk_table$end),
                  sep="\t"),
            sep="\n",file = separator,append=F)
        
      #--------- create Highlight.conf
      setwd("../conf/")
      highl <- paste(SRA,"-",nam, "_highlights.conf", sep = "")
      
      cat(paste("<highlights>",sep=""),sep="\n",file = highl,append=TRUE)
      
      cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("file = data/",cover1, sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r0   = 0.47r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r1   = 0.83r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      
      cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("file = data/",gene.blocks , sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r0   = 0.50r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r1   = conf(.,r0)+0.03r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      
      cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("file = data/",highlight1 , sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r0   = 0.83r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r1   = conf(.,r0)+0.03r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      
      cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("file = data/",separator , sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r0   = 0.89r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r1   = conf(.,r0)+0.005r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("fill_color = vlgrey", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      
      cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("file = data/",separator , sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r0   = 0.79r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("r1   = conf(.,r0)+0.005r", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("fill_color = vlgrey", sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
      cat(paste("</highlights>",sep=""),sep="\n",file = highl,append=TRUE)
      
      #--------- create ideogram.conf
      ideo <- paste(SRA,"-",nam, "_ideogram.conf", sep = "")
      
      cat(paste("<ideogram>",sep=""),sep="\n",file = ideo,append=TRUE)
      
      cat(paste("<spacing>",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("default = 0u",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("break   = 0u",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("</spacing>",sep=""),sep="\n",file = ideo,append=TRUE)
      
      cat(paste("thickness        = 20p",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("stroke_thickness = 2",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("stroke_color     = black",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("fill             = yes",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("fill_color       = black",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("radius         = 0.85r",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("show_label     = no",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("label_font     = default",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("label_radius   = dims(ideogram,radius) + 0.05r",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("label_size     = 36",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("label_parallel = yes",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("label_case     = upper",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("band_stroke_thickness = 2",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("show_bands            = yes",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("fill_bands            = yes",sep=""),sep="\n",file = ideo,append=TRUE)
      cat(paste("</ideogram>",sep=""),sep="\n",file = ideo,append=TRUE)
      
      #--------- create ticks_plasmid.conf
      ticks <- paste(SRA,"-",nam, "_ticks.conf", sep = "")
      
      cat(paste("show_ticks          = yes",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("show_tick_labels    = yes",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("show_grid          = no",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_start         = dims(ideogram,radius_inner)-0.5r",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_end           = dims(ideogram,radius_inner)",sep=""),sep="\n",file = ticks,append=TRUE)
      
      cat(paste("<ticks>",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("skip_first_label     = no",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("skip_last_label      = no",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("radius               = dims(ideogram,radius_outer)",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("tick_separation      = 2p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("min_label_distance_to_edge = 0p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("label_separation = 5p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("label_offset     = 5p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("label_size = 8p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("multiplier = 1",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("color = black",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("thickness = 3p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("size      = 20p",sep=""),sep="\n",file = ticks,append=TRUE)
      
      cat(paste("<tick>",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("size           = 10p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("spacing        = 0.5u",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("color          = black",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("show_label     = yes",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("label_size     = 30p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("format         = %s",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid           = no",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_color     = lgrey",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_thickness = 1p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("</tick>",sep=""),sep="\n",file = ticks,append=TRUE)
      
      cat(paste("<tick>",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("size           = 15p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("spacing        = 5u",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("color          = black",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("show_label     = yes",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("label_size     = 45p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("format         = %s",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid           = yes",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_color     = lgrey",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("grid_thickness = 1p",sep=""),sep="\n",file = ticks,append=TRUE)
      cat(paste("</tick>",sep=""),sep="\n",file = ticks,append=TRUE)
      
      cat(paste("</ticks>",sep=""),sep="\n",file = ticks,append=TRUE)
      
      #--------- create Final.conf
      final <- paste(SRA,"-",nam, ".conf", sep = "")
      
      cat(paste("<<include etc/colors_fonts_patterns.conf>>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<<include ", ideo, ">>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<<include ", ticks, ">>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("karyotype = data/",karyo,sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<image>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<<include etc/image.conf>>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("</image>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<<include ",highl,">>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("chromosomes_units           = 1000",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("chromosomes_display_default = yes",sep=""),sep="\n",file = final,append=TRUE)
      
      cat(paste("<plots>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("<plot>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("type  = text",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("file  = data/",label.locus,sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("color = black",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("r1    = 0.95r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("r0    = 0.55r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_size = 26",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_font = light",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("padding    = 5p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("rpadding   = 5p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("show_links     = yes",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("link_dims      = 5p,4p,8p,4p,0p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("link_thickness = 1p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("link_color     = dgrey",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_snuggle        = yes",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("max_snuggle_distance = 2r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_sampling     = 1",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_tolerance    = 0.25r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_link_overlap_test      = yes",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_link_overlap_tolerance = 2p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_refine                 = no",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("</plot>",sep=""),sep="\n",file = final,append=TRUE)
      
      cat(paste("<plot>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("type  = text",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("file  = data/",label.genes,sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("color = black",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("r1    = 0.45r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("r0    = 0.30r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_size = 26",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_font = bold",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("padding    = 5p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("rpadding   = 5p",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("show_links     = no",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("label_snuggle        = yes",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("max_snuggle_distance = 2r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_sampling     = 1",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_tolerance    = 0.25r",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("snuggle_refine                 = no",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("</plot>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("</plots>",sep=""),sep="\n",file = final,append=TRUE)
      
      cat(paste("<<include etc/housekeeping.conf>>",sep=""),sep="\n",file = final,append=TRUE)
      cat(paste("data_out_of_range* = trim",sep=""),sep="\n",file = final,append=TRUE)
      
      setwd("../../../../")

}
}

#--------- create  Circos.sh for perl using local:lib 
setwd("circos/circos-0.69-6/circos_plasmid/conf/")
my_save_sh <- list.files(pattern = "\\_highlights.conf")
my_save_sh <- sapply(my_save_sh, function(x) gsub(".*-","",x))
my_save_sh <- sapply(my_save_sh, function(x) gsub("_highlights.conf","",x))


Script_locallib  = "perl_local_circos.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script_locallib,append=TRUE)

for (sr in 1:length(SRA_names)) {
    for (i in 1:length(my_save_sh)) {
      cat(paste("perl -Mlocal::lib ../../bin/circos -conf ", paste(SRA_names[sr],"_",my_save_sh[i], ".conf", sep = ""), " -outputfile ", paste(SRA_names[sr],"_",my_save_sh[i], sep= ""), sep = ""),sep="\n",file=Script_locallib,append=TRUE)  
    }
}

setwd("../../../../")
