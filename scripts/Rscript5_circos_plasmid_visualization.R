
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/")

### Inputs ###
#============

output_list #from Rscripts2_magicBlast_output

#============


### Code ###
#============


#--------- load gff
setwd("gff3_plasmid_files/")
my_gff_files <- list.files(pattern = "\\_Gff.csv")
my_gff <- sapply(my_gff_files, function(x) gsub("\\..*","",x))

x <- which(names(output_list) %in% my_gff)
save <- names(output_list[x])

for (r in 1:length(x)) {
  
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/gff3_plasmid_files/")
nam <- names(output_list[x[r]])
col.start <- which(colnames(output_list[[x[r]]]) == "contigs_number")
col.end <- ncol(output_list[[x[r]]])
new.df <- output_list[[x[r]]][, col.start:col.end]

nodupl.df <- new.df[!duplicated(new.df[,1]),]

my_gff_table <- read.csv(my_gff_files[[r]])
size.plas <- max(my_gff_table$end)



#--------- create Karyotype
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/circos/circos-0.69-6/staph_plasmids/data/staph/")
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
for (i in 1:nrow(my_gff_table)) {
  
  cat(paste("chr1", my_gff_table$start[i], 
            my_gff_table$end[i],
            paste("fill_color=chr", sample(1:20,1), sep=""),
            sep="\t"),
      sep="\n",file = gene.blocks,append=T)
}
temp <- read.delim(paste(nam, "_gene_block.txt", sep = ""), sep = "\t", header = F)
temp <- temp[-1,]
write.table(temp, paste(nam, "_gene_block.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")


#--------- create label_locus.txt
label.locus <- paste(nam, "_locus_labels.txt", sep = "")
for (i in 1:nrow(my_gff_table)) {
  
  cat(paste("chr1", 
            my_gff_table$start[i], 
            my_gff_table$end[i],
            my_gff_table$GenBank_locus_tag[i],
            sep="\t"),
      sep="\n",file = label.locus,append=T)
}
temp <- read.delim(paste(nam, "_locus_labels.txt", sep = ""), sep = "\t", header = F)
temp <- temp[-1,]
write.table(temp, paste(nam, "_locus_labels.txt", sep = ""), quote = F, row.names = F, col.names = F, sep = "\t")


#--------- create label_genes.txt
label.genes <- paste(nam, "_gene_labels.txt", sep = "")
for (i in 1:nrow(my_gff_table)) {
  
  cat(paste("chr1", 
            my_gff_table$start[i], 
            my_gff_table$end[i],
            my_gff_table$Gene_name[i],
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
            max(my_gff_table$end),
            sep="\t"),
      sep="\n",file = separator,append=F)
  

#--------- create Highlight.conf
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/circos/circos-0.69-6/staph_plasmids/conf/")
highl <- paste(nam, "_highlights.conf", sep = "")

cat(paste("<highlights>",sep=""),sep="\n",file = highl,append=TRUE)

cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("file = data/staph/",cover1, sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r0   = 0.47r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r1   = 0.83r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)

cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("file = data/staph/",gene.blocks , sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r0   = 0.50r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r1   = conf(.,r0)+0.03r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)

cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("file = data/staph/",highlight1 , sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r0   = 0.83r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r1   = conf(.,r0)+0.03r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)

cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("file = data/staph/",separator , sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r0   = 0.89r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r1   = conf(.,r0)+0.005r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("fill_color = vlgrey", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)

cat(paste("<highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("file = data/staph/",separator , sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r0   = 0.79r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("r1   = conf(.,r0)+0.005r", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("fill_color = vlgrey", sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlight>",sep=""),sep="\n",file = highl,append=TRUE)
cat(paste("</highlights>",sep=""),sep="\n",file = highl,append=TRUE)


#--------- create Final.conf
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/circos/circos-0.69-6/staph_plasmids/conf/")
final <- paste(nam, ".conf", sep = "")

cat(paste("<<include etc/colors_fonts_patterns.conf>>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<<include ideogram.conf>>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<<include ticks_plasmid.conf>>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("karyotype = data/staph/",karyo,sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<image>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<<include etc/image.conf>>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("</image>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<<include ",highl,">>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("chromosomes_units           = 1000",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("chromosomes_display_default = yes",sep=""),sep="\n",file = final,append=TRUE)

cat(paste("<plots>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("<plot>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("type  = text",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("file  = data/staph/",label.locus,sep=""),sep="\n",file = final,append=TRUE)
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
cat(paste("snuggle_refine                 = yes",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("</plot>",sep=""),sep="\n",file = final,append=TRUE)

cat(paste("<plot>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("type  = text",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("file  = data/staph/",label.genes,sep=""),sep="\n",file = final,append=TRUE)
cat(paste("color = black",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("r1    = 0.45r",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("r0    = 0.35r",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("label_size = 26",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("label_font = bold",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("padding    = 5p",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("rpadding   = 5p",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("show_links     = no",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("label_snuggle        = yes",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("max_snuggle_distance = 2r",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("snuggle_sampling     = 1",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("snuggle_tolerance    = 0.25r",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("snuggle_refine                 = yes",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("</plot>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("</plots>",sep=""),sep="\n",file = final,append=TRUE)

cat(paste("<<include etc/housekeeping.conf>>",sep=""),sep="\n",file = final,append=TRUE)
cat(paste("data_out_of_range* = trim",sep=""),sep="\n",file = final,append=TRUE)

}

#--------- create  Circos.sh 
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/scripts/")
Script = "circos.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
for (i in 1:length(x)) {
  cat(paste("perl ../../bin/circos -conf ", paste(save[i], ".conf", sep = ""), " -outputfile ", save[i], sep = ""),sep="\n",file=Script,append=TRUE)  
}
# chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid/scripts/circos.sh


#--------- create  Circos.sh for perl using local:lib
setwd("C:/Users/copinrj/Desktop/NCBI/Got_plasmid_project/scripts/")
Script_locallib  = "perl_local_circos.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script_locallib,append=TRUE)
for (i in 1:length(x)) {
  cat(paste("perl -Mlocal::lib ../../bin/circos -conf ", paste(save[i], ".conf", sep = ""), " -outputfile ", save[i], sep = ""),sep="\n",file=Script_locallib,append=TRUE)  
}
# chmod +x /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/got_plasmid/scripts/perl_local_circos.sh

