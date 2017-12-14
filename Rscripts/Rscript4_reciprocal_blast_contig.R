
### Functions ###
#===============
readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n")
  chars = substr(row,1,1)
  base = chars==">"
  seq = row[base]
  return(seq)
}
#===============

### Inputs ###
#===============
setwd("plasmids/contigs/")
contig_files = list.files(pattern = "\\_matching_contigs.fasta")
contig_names = sapply(contig_files, 
                      function(s) {s = gsub("_matching_contigs.fasta", "", s)})
SRA <- "SRR6227128"

#===============

### Code ###
#===========
#--------- identify empty blast results
setwd("reciprocal_blast/")

empty <- list()
my_blast_files <- list.files(pattern = "\\vs")
for (h in 1:length(my_blast_files)) {
  info = file.info(my_blast_files[h])
  empty[[h]] = rownames(info[info$size == 0, ])
}
names(empty) <- my_blast_files
empty <- empty[lapply(empty, length)>0]
names(empty)

no.recipro <- as.data.frame(names(empty))
colnames(no.recipro) <- "file"
no.recipro$genome1 <- no.recipro$file
no.recipro$genome2 <- no.recipro$file
no.recipro$genome1 <- sapply(no.recipro$genome1, function(s) gsub("vs.*","",s))
no.recipro$genome2 <- sapply(no.recipro$genome2, function(s) gsub(".*vs","",s))
no.recipro$genome2 <- sapply(no.recipro$genome2, function(s) gsub(".txt","",s))

#--------- work with the rest

recipro_names <- my_blast_files[- which(my_blast_files %in% names(empty))]
my_recipro_list <- lapply(recipro_names, function(i){fread( i, sep = "\t", header=F, data.table = F, fill = T)})
names(my_recipro_list) <- recipro_names



for (i in 1:length(contig_names)) {
  for (j in 1:length(contig_names)) {
    
  
    blast.file.name1 = paste(contig_names[j], "vs", contig_names[i], ".txt", sep = "")
    blast.file.name2 = paste(contig_names[i], "vs", contig_names[j], ".txt", sep = "")
    

    
    if(blast.file.name1 %in% recipro_names ==T){if(blast.file.name2 %in% recipro_names==T){
    Blast1 <- my_recipro_list[[blast.file.name1]]
    colnames(Blast1) <- c("Query_id_a", "Subject_id_b", "%_identity_a", "alignment_length_a", "mismatches_a", "gap_openings_a", "q_start_a", "q_end_a", "s_start_a", "s_end_a", "e_value_a", "bit_score_a", "Subject_b_GI", "Subject_b_Title" )
    
    # select alignments with length of at least 50 nu
    
    Blast1.align50nu <- Blast1[which(Blast1$alignment_length_a >= 50),c(1,2,3,4,11,12,14)]
    
    # for each query, extract the blast output for their corresponding target/s 
    # sort ascending based on bit_score
    # store matrices in a list
    
    unique.row <- unique(Blast1.align50nu$Query_id_a)
    atBeginning <- Sys.time()
    then <- Sys.time()
    list.best.hits1 <- list()
    list.best.hits2 <- list()
    
    for (k in 1:length(unique.row)){
      d <- Blast1.align50nu[which(Blast1.align50nu$Query_id_a == unique.row[k]),]
      list.best.hits1[[k]] <- d[order(d$bit_score_a, decreasing = T),]
      list.best.hits2[[k]] <-  list.best.hits1[[k]][1:4,]
      names(list.best.hits1)[k] <- paste(list.best.hits2[k])
    }
    print(Sys.time()-then)
    
    # for each query protein, extract the Best Hit of the Query (BHQ) i.e. target protein/s with the best bit_score  (store as list)
    
    atBeginning <- Sys.time()
    then <- Sys.time()
    targets.a <- list()
    
    for (tg in names(list.best.hits1)){
      targets.a[[tg]] <- list.best.hits1[[tg]][which(list.best.hits1[[tg]]$bit_score_a == max(list.best.hits1[[tg]]$bit_score_a)),2]
    }
    print(Sys.time()-then)
    
  
    
    ######## READ in the blast output file for the genome2 #######
    
    Blast2 <- my_recipro_list[[blast.file.name2]]
    colnames(Blast2) <- c("Query_id_b", "Subject_id_a", "%_identity_b", "alignment_length_b", "mismatches_b", "gap_openings_b", "q_start_b", "q_end_b", "s_start_b", "s_end_b", "e_value_b", "bit_score_b", "Subject_a_GI", "Subject_a_Title")
    
    # select alignments with length of at least 50 nu
    
    Blast2.align50nu <- Blast2[which(Blast2$alignment_length_b >= 50),c(1,2,3,4,11,12,14)]
    Blast2.align50nu <- Blast2.align50nu[order(Blast2.align50nu$bit_score_b, decreasing = T),]
    
    # for each query, extract the blast output for their corresponding BHQ
    # sort decreasing based on bit_score
    # store matrices in a list
    
    atBeginning <- Sys.time()
    then <- Sys.time()
    targets.b <- list()
    
    for (q in names(targets.a)){
      d <- Blast2.align50nu[which(Blast2.align50nu$Query_id_b %in% targets.a[[q]]),]
      targets.b[[q]] <- d[order(d$bit_score_b, decreasing = T),]
    }
    print(Sys.time()-then)
    
    # for the BHQ of each query protein, if the BHT (protein/s with the smallest E-value) match the query protein => assign orthology as RBH
    # store matrices in a list
    
    atBeginning <- Sys.time()
    then <- Sys.time()
    final.list <- list()
    
    for (r in names(targets.b)){
      
      
      
      maxi <- targets.b[[r]][which(targets.b[[r]]$bit_score_b == max(targets.b[[r]]$bit_score_b)),]
      sub.maxi <- targets.b[[r]][-which(targets.b[[r]]$bit_score_b == max(targets.b[[r]]$bit_score_b)),]
      sub.maxi <- sub.maxi[which(sub.maxi$bit_score_b == max(sub.maxi$bit_score_b)),]
      maxi <- rbind(maxi, sub.maxi)
      final.list[[r]] <- maxi[which(maxi$Subject_id_a %in% unique(list.best.hits1[[paste(r)]]$Query_id_a)),]

    }
    print(Sys.time()-then)
    
    # store the the RBH as matrix 
    
    atBeginning <- Sys.time()
    then <- Sys.time()
    final.final.list <- do.call("rbind",final.list)
    rownames(final.final.list) <- NULL
    print(Sys.time()-then)
    
    # Time difference of 3.460173 mins
    

    
    final.final.list$genome_query<- gsub(".txt","",gsub(".*vs","",blast.file.name1))
    final.final.list$genome_subject  <- gsub("vs.*","",blast.file.name1)
    final.final.list$genome_file <- blast.file.name2


    # save the RBH matrix
    
    
    write.table(final.final.list, paste(contig_names[j], "___",SRA,"___", contig_names[i], ".txt", sep=""), quote=FALSE, sep = "\t", row.names = F)
    }}
  }
}


files = list.files(pattern = SRA)
reciprocal.df <- lapply(files, function(i){fread( i, sep = "\t", header=T, data.table = F, fill = T)})

reciprocal.df <- do.call(rbind,reciprocal.df)
reciprocal.df$contig.query <- paste(reciprocal.df$genome_query, reciprocal.df$Query_id_b)
reciprocal.df$contig.subject <- paste(reciprocal.df$genome_subject, reciprocal.df$Subject_id_a)
unique.contig <- unique(reciprocal.df$contig.query)
unique.contig <- unique.contig[order(unique.contig)]

cross.table <-  matrix(NA, length(unique(reciprocal.df$contig.query)),length(unique(reciprocal.df$contig.query)))
cross.table <- as.data.frame(cross.table)
row.names(cross.table) <- unique(reciprocal.df$contig.query[order(reciprocal.df$contig.query)])
colnames(cross.table) <-  unique(reciprocal.df$contig.query[order(reciprocal.df$contig.query)])


for (i in 1:nrow(cross.table)) {
  
  f <- reciprocal.df$contig.subject[which(reciprocal.df$contig.query %in% unique.contig[i])]
  cross.table[i,which(colnames(cross.table) %in% f)] <- 1
  cross.table[which(rownames(cross.table) %in% f),i] <- 1
}

for(i in 1:ncol(cross.table)) {
  cross.table[which(is.na(cross.table[,i]) ==T),i] <- 0 
}

write.csv(cross.table, paste("../../../outputs/", SRA, "_contig_cross_table.csv", sep = ""))

x <- as.data.frame(colSums(cross.table))


colnames(x) <- "Number of matching contigs"
x$plasmid <- sapply(row.names(x), function(s)gsub(" .*", "", s))
x <- x[c(2,1)]

y <- as.data.frame(matrix(NA, nrow(x), max(x$`Number of matching contigs`)))
y$plasmid <- x$plasmid
y$contig <- rownames(x)
y$`Number of matching contigs` <- x$`Number of matching contigs`
y$`Number of matching plasmid` <- NA
y <- y[c((ncol(y)-3):ncol(y), 1:(ncol(y)-4))]

for(i in 1: nrow(y)){
  
  f <- which(colnames(cross.table) == y$contig[i])
  f1 <- row.names(cross.table)[which(cross.table[,f] == 1)]
  f2 <- sapply(f1, function(x) gsub(" .*", "",x))
  f3 <- unique(f2)
  y$`Number of matching plasmid`[i] <- length(f3)
  y[i,5:(length(f3)+4)] <- f3
  
}

for(i in 1:ncol(y)) {
  y[which(is.na(y[,i]) ==T),i] <- "-" 
}

write.csv(y, paste("../../../outputs/", SRA, "_summary_table.csv", sep = ""), row.names = F)

setwd("../../../")
