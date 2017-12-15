



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
SRA <- scan("SRA/SRA_ID.txt", what=character(0),sep="\n")
setwd("plasmids/magic_output/")

for (sr in 1:length(SRA)) {
  
    tg.SRA = SRA[sr]
    
    my_magic_files <- list.files(pattern = tg.SRA)
    my_magic_names <- gsub("_magic.txt","",my_magic_files)
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
        my_plas_name <- plasmid_names[o]
        my_magic_name <- my_magic_names[o]
        ref <- readFastaRef(paste("../fasta/", my_plas_name, ".fasta", sep = ""))
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
        names(output_list)[o] <- plasmid_names[o]
    
        #--------- create contig.fasta file
        outfile <- paste("../contigs/", my_magic_name,"_matching_contigs.fasta", sep = "")
        
        for(i in 1: nrow(coord)){
            a <- ref[c(coord$start[i]:coord$end[i])]
            b <- paste(a,collapse="")
            cat(paste(">contig_",i, sep=""),sep="\n",file = outfile,append=TRUE)
            cat(paste(b, sep=""),sep="\n",file = outfile,append=TRUE)
        }
        
    }
    
    #--------- Create file with contig sizes for each plasmid
    contig_sizes <- list()
    
        for(d in 1: length(output_list)){
        u.contig <- as.data.frame(unique(output_list[[d]]$contigs_number))
        vect <- c(NA)
    
            for(x in 1:nrow(u.contig)){
            size <- max(output_list[[d]]$overl_size[which(output_list[[d]]$contigs_number %in% as.character(u.contig[x,1]))])
            vect <- append(vect, size)
            }
        
        vect <- vect[-1]
        u.contig$size <- vect
        u.contig$genome <- names(output_list)[[d]]
        colnames(u.contig) <- c("contig", "size", "genome")
        contig_sizes[[d]] <- u.contig
    }
    
    size.contig <- do.call(rbind, contig_sizes)
    
    
    
    #--------- Create a centrilized table to identify full plasmid signature 
    m <- as.data.frame(matrix(NA,length(output_list),5))
    colnames(m) <- c("plasmid_ID","size_plasmid","Number_of_contigs","size_overlap","Overlapping_portion(%)")
    
    for (r in 1:length(output_list)){
        m$plasmid_ID[r] <-  names(output_list)[r]
        m$Number_of_contigs[r] <- max(output_list[[r]]$overlapping_contigs)
        m$size_plasmid[r] <- max(output_list[[r]]$ref_size)
        m$size_overlap[r] <- sum(unique(output_list[[r]]$overl_size))
        m$`Overlapping_portion(%)`[r] <- max(output_list[[r]]$overall_coverage)
    }
    
    #--------- Save the table
    write.csv(m,paste("../../outputs/",SRA[sr],"_got_plamid.csv", sep = ""), row.names = F)
}
setwd("../../")
#==========================================================================================================================================




