
### Inputs ###
#============

output_list #from Rscripts2_magicBlast_output

#============


### Code ###
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
      karyo <- paste(nam, "_karyotpye.txt", sep = "")
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
      highlight1 <- paste(nam, "_highlight1.txt", sep = "")
            for (i in 1:round(nrow(nodupl.df))) {
                
                cat(paste("chr1", nodupl.df$overl_coord_start[i], 
                          nodupl.df$overl_coord_end[i],
                          paste("fill_color=chr", sample(1:20,1), sep=""), 
                          sep="\t"),
                    sep="\n",file = highlight1,append=T)
            }
      
      #--------- create cover.txt
      cover1 <- paste(nam, "_cover1.txt", sep = "")
            for (i in 1:round(nrow(nodupl.df))) {
              
              cat(paste("chr1", nodupl.df$overl_coord_start[i], 
                        nodupl.df$overl_coord_end[i],
                        paste("fill_color=vlgrey",sep=""), 
                        sep="\t"),
                  sep="\n",file = cover1,append=T)
            }
      
      #--------- create gene_block.txt
      gene.blocks <- paste(nam, "_gene_block.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        paste("fill_color=chr", sample(1:20,1), sep=""),
                        sep="\t"),
                  sep="\n",file = gene.blocks,append=T)
            }
      temp <- read.delim(paste(nam, "_gene_block.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      write.table(temp, paste(nam, "_gene_block.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create label_locus.txt
      label.locus <- paste(nam, "_locus_labels.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", 
                        my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        my_gbk_table$locus_tag[i],
                        sep="\t"),
                  sep="\n",file = label.locus,append=T)
            }
      temp <- read.delim(paste(nam, "_locus_labels.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      write.table(temp, paste(nam, "_locus_labels.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create label_genes.txt
      label.genes <- paste(nam, "_gene_labels.txt", sep = "")
            for (i in 1:nrow(my_gbk_table)) {
              
              cat(paste("chr1", 
                        my_gbk_table$start[i], 
                        my_gbk_table$end[i],
                        my_gbk_table$gene[i],
                        sep="\t"),
                  sep="\n",file = label.genes,append=T)
            }
      temp <- read.delim(paste(nam, "_gene_labels.txt", sep = ""), sep = "\t", header = F)
      temp <- temp[-1,]
      s <- which(is.na(temp$V4) == T)
            if (length(s)>0) {
              temp <- temp[-which(is.na(temp$V4) == T),]
            }
      write.table(temp, paste(nam, "_gene_labels.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")
      
      #--------- create separator.txt
      separator <- paste(nam, "_separator.txt", sep = "")
        cat(paste("chr1", 
                  0, 
                  max(my_gbk_table$end),
                  sep="\t"),
            sep="\n",file = separator,append=F)
        
      #--------- create Highlight.conf
      setwd("../conf/")
      highl <- paste(nam, "_highlights.conf", sep = "")
      
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
      ideo <- paste(nam, "_ideogram.conf", sep = "")
      
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
      ticks <- paste(nam, "_ticks.conf", sep = "")
      
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
      final <- paste(nam, ".conf", sep = "")
      
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


#--------- create  Circos.sh for perl using local:lib 
setwd("circos/circos-0.69-6/circos_plasmid/conf/")
my_save_sh <- list.files(pattern = "\\_highlights.conf")
my_save_sh <- sapply(my_save_sh, function(x) gsub("_highlights.conf","",x))


Script_locallib  = "perl_local_circos_temp.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script_locallib,append=TRUE)
for (i in 1:length(my_save_sh)) {
  cat(paste("perl -Mlocal::lib ../../bin/circos -conf ", paste(my_save_sh[i], ".conf", sep = ""), " -outputfile ", my_save_sh[i], sep = ""),sep="\n",file=Script_locallib,append=TRUE)  
}
# chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid/scripts/perl_local_circos.sh
setwd("../../../../")
