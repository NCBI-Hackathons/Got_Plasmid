# Got Plasmids?
## (or phages, repetitive elements or polintons?)

## Goal
Defining the clinically relevant plasmids of Staphylococcus aureus and other pathogenic bacteria.

## Challenge: 
Identify plasmids using Next generation sequencing  and SRA data.
Mobile genetic elements, including plasmids and phages, are ambiguous as rich in repeats and low complexity sequences.

AIMS:
- Identification of unique plasmids in individual patient samples
- Differentiation between screening (nose swabs) and lesions
- Building a rapid pipeline to differentiate between plasmids in these differential sites.  
- Extending that to other bacterial species

 ## Visualization using Circos
 
 Contigs from same SRA dataset mapped against 3 different S. aureus plasmids.
 Contigs are represented in the first colored circle.
 Genes are in the second colored circled.
 Genes names and locus tag are indicated when available.
 
 Plasmid present in SRA - 2 contigs (blue, pink):
 ![My image](https://github.com/NCBI-Hackathons/Pathogenic_Pangenomes/blob/master/images/AP003139.png)
 
 Plasmid partially present:
 ![My image](https://github.com/NCBI-Hackathons/Pathogenic_Pangenomes/blob/master/images/GQ900377.1.png)
 
 Plasmid absent in SRA - no contig homology:
 ![My image](https://github.com/NCBI-Hackathons/Pathogenic_Pangenomes/blob/master/images/GQ900379.1.png)
 

## Dependencies ###

magicBlast

Blast

Perl

Circos

R

## Setup ###

- git clone https://github.com/NCBI-Hackathons/Got_plasmid.git

###CIRCOS###
- Circos runs in perl
- dowload circos-0.69-6.tgz from http://circos.ca/software/download/circos/
- mv circos-0.69-6.tgz Path_to/Got_plasmid/circos/
- cd Path_to/Got_plasmid/circos/
- tar -xzvf circos-0.69-6.tgz
- rm *.tgz

- mkdir circos-0.69-6/circos_plasmid/
- mkdir circos-0.69-6/circos_plasmid/conf/
- mkdir circos-0.69-6/circos_plasmid/data/

- check perl modules: perl -Mlocal::lib circos-0.69-6/bin/circos -modules
- if missing modules
- On Unix, you can install perl modules locally : 
- perl -MCPAN -Mlocal::lib -e shell
- cpan[1]>install module_name (example cpan[1]>install Math::Bezier)
- For other machines, read instructions: http://circos.ca/documentation/tutorials/configuration/perl_and_modules/

###R### 
- cd Path_to/Got_plasmid/
- wget http://cran.rstudio.com/src/base/R-3/R-3.4.1.tar.gz 
- tar xvf R-3.4.1.tar.gz $ cd R-3.4.1 $ ./configure --prefix=$HOME/R 
- make && make install 
- export PATH=$PATH:/HOME/R/bin

run R
- R
- source("https://bioconductor.org/biocLite.R")
- biocLite("GenomicRanges")
- biocLite("IRanges")
- biocLite("Biostrings")
- install.packages("data.table")
- install.packages("reutils")
- install.packages("devtools")
- install.packages("biofiles")
- install.packages("Biostrings")
- devtools::install_github("gschofl/biofiles")


## WorkFlow

 1.  Extract all genomic gapless genome fasta and gff.
 2.  Remove plasmid sequences from genomes.
 3.  Retreive plasmid fasta with eUtils.
 4.  Create customized blast databases.
 5.  Use magicBlast on one SRRA versus individual plasmid databases.
 6.  Create individual contig.fasta files and generate a table with the % of plasmid sequences covered by contigs.
 7.  Download and parse plasmid GenBank genome files.
 8.  Generate plasmid and contig visualization using Circos.
 9.  Create customized contig db and blast contigs reciprocally.
 10. Identify contigs matching each other and other plasmids.
 11. BLAST contigs against gapless genome databases.
 12. Identify contigs matching gapless genome sequences.


## STEPS 


### Step 1.

      # Extract all genomic gapless genome fasta and gff.

      # from NCBI website, go to Assembly database:
      # staphylococcus aureus[Organism] 
      # Filters: Status = Latest
      # Filters: Assembly level = complete genomes
      # download Assembly IDs
      # extract first column and remove column names
      # save as "gapless_genome_assemblies.txt" in Path_to/Got_plasmid/gapless_genomes/assembly_ID/

      ## ***command line***
      ## cd Path_to/Got_plasmid/gapless_genomes/assembly_ID/
      ## bash genome_getter.sh
      #
      # genome_getter.sh:
      # -----------------
      # #!/bin/bash
      # for i in `cat gapless_genome_assemblies.txt`; do wget `esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`; done
      # gunzip *.gz
      # cat *.fna > gapless_genome_plasmid.fasta
      # mv gapless_genome_plasmid.fasta ../fasta/.

      ## ***command line***
      ## cd Path_to/Got_plasmid/gapless_genomes/assembly_ID/
      ## bash gff_getter.sh
      #
      # gff_getter.sh:
      # -------------
      # #!/bin/bash
      # for i in $(cat gapless_genome_assemblies.txt); do wget $(esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.gff.gz"}'); done
      # gunzip *.gz
      # mv *.gff ../gff/.


# Step 2. 
# ------------------
      # Remove plasmid sequences from genomes.

      ## ***command line***
      ## cd Path_to/Got_plasmid/

      ### ./R
      ### Use Rscript1_split_genome_from_plasmid_sequences.R 
      ### q()

# Step 3. 
# ------------------
      # Retreive plasmid fasta with eUtils.

      # from NCBI website, go to Nucleotide database:
      # plasmid[title] AND staphylococcus[title]
      # Filters: Species = bacteria
      # Filters: Molecule types = genomic RNA/DNA
      # Filters: Genetic compartments = Plasmid
      # download accession table
      # save as "accession_plasmids.txt" in /plasmids/assembly_ID/

      ## ***command line***
      ## cd Path_to/Got_plasmid/

      ### ./R
      ### use esearch_plasmid.R
      ### q()

      ## ***command line***
      ## cd Path_to/Got_plasmid/plasmids/assembly_ID/
      ## bash esearch_plasmid.sh

# Step 4. 
# ------------------
      # Create customized blast databases.

      ## ***command line***
      ## cd Path_to/Got_plasmid/gapless_genomes/fasta/
      ## bash makeblastdb_gapless_genomes.sh
      ## cd Path_to/Got_plasmid/

      ### ./R
      ### use makeblastdb_individual_plasmid.R
      ### q()

      ## ***command line***
      ## cd PAth_to/Got_plasmid/plasmids/fasta/
      ## bash makeblastdb_individual_plasmid.sh

      
# Step 5. 
# ------------------
      # Use magicBlast on one SRRA versus individual plasmid databases.

      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use magicblast_plasmid.R
      ### q()

      ## ***command line***
      ## cd PAth_to/Got_plasmid/plasmids/magic_output/
      ## bash magicblast_plasmid.sh


# Step 6. 
# ------------------
      # Create individual contig.fasta files and generate a table with the % of plasmid sequences covered by contigs.

      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use Rscript2_magicBlast_output.R
      ### q()


# Step 7. 
# ------------------
      # Download and parse plasmid GenBank genome files.
      
      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use Rscript3_GenBank_parser.R
      ### q()


# Step 8. 
# ------------------
      # Generate plasmid and contig visualization using Circos.
      
      ## ***command line***
      ## cd Path_to/Got_plasmid/

      ### ./R
      ### use Rscript4_circos_plasmid_visualization
      ### q()


      ## ***command line***
      ## cd PAth_to/Got_plasmid/circos/circos-0.69-6/circos_plasmid/conf/
      ## bash perl_local_circos_temp.sh
      ## mv *.png ../../../../outputs/.


# Step 9. 
# ------------------
      # Create customized contig db and blast contigs reciprocally.

      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use contigs.R
      ### q()

      ## ***command line***
      ## cd PAth_to/Got_plasmid/plasmids/contigs/
      ## bash Formatblast_contigs.sh
      ## bash reciprocal_contig_Blast.sh
      
      
# Step 10. 
# ------------------
      # Identify contigs matching each other and other plasmids.

      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use Rscript5_reciprocal_blast_contig.R
      ### q()


# Step 11. 
# ------------------
      # BLAST contigs against gapless genome databases.
      
      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use blast_contigs_vs_gapless_genomes.R
      ### q()
      
      ## ***command line***
      ## cd PAth_to/Got_plasmid/plasmids/contigs/
      ## ./blast_contigs_vs_gapless_genomes.sh



# Step 12. 
# ------------------
      # Identify contigs matching gapless genome sequences.
    
      ## ***command line***
      ## cd Path_to/Got_plasmid/
      
      ### ./R
      ### use Rscript6_contig_match_in_genomes.R
      ### q()


# ==============


## Backgound and Significance
### Staphylococcus aureus genomes and mobile genetic elements (MGE).
The proliferation of S. aureus genomic sequences in public databases reflects strong interest in understanding S. aureus genome diversity and evolution. With an average of 2,800 coding sequences, it is estimated that 44% of S. aureus genes are NOT shared by all S. aureus strains. These genes constitute the  ‘accessory genome’, which is variable between strains and mostly made of mobile genetic elements (MGE) enriched in hypothetical proteins. 

Staphylococcal MGE encompass any intra- or extra-chromosomal DNA segment that can be independently mobilized within or between S. aureus cells. It includes plasmids, transposons, integrons, genomic islands, S. aureus pathogenicity islands (SaPIs), integrative conjugative elements, staphylococcal chromosome cassettes, and phages. Together, phages and plasmids are the main source of MGE diversity among S. aureus strains. 

### Why study Staphylococcus aureus plasmid diversity?
MGE discovery and characterization are important goals for clinical genomic analysis because almost all S. aureus strains harbor at least one plasmid with potentially syndrome- and tissue-specific functions. 
 
### Plasmid diversity in NCBI.
 As of August 2017,  327 unique plasmid sequences have been identified and deposited in the US National Center for Biotechnology and Information (NCBI) database.  With the ease and speed of whole genome sequencing, new MGE are discovered every day, highlighting the impressive breadth of S. aureus plasmid diversity. There are ~40,000 primary (unannotated) S. aureus datasets in the Sequence Read Archive (SRA). 

### The challenge to study MGE diversity. 
The extent and importance of plasmid contribution to S. aureus pathogenesis is largely under-appreciated. This is mainly due to the complications inherent to their identification, analysis and characterization. 
 
Plasmid are rich in repetitive sequences. High level of sequence identity facilitates genetic recombination and contributes to the emergence of mosaic MGE. As such, MGE are ever-changing and can be hard to identify. Moreover, nomenclature in public databases is constantly evolving and inconsistency in annotation among MGE is common and complicates functional inter- and intra-species analyses. 


