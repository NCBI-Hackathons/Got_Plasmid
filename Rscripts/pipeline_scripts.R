#
# -- Files
# -- Shell Sripts
#    | genome_getter.sh
#    | gff_getter.sh
#    | esearch_plasmid.sh
#    | makeblastdb_individual_plasmid.sh
#    | makeblastdb_all_plasmids.sh
#    | makeblastdb_gapless_genomes.sh
#    | magicblast_plasmid.sh
#    | BLAST_contig_gainst_all_plasmids.sh
#    | BLAST_contig_against_gapless_genomes.sh
# -- DeBugging/Notes
#
##  ##

### WorkFlow ###
# ===========

# 1.  Extract all genomic gapless genome fasta and gff.
# 2.  Remove plasmid sequences from genomes.
# 3.  Retreive plasmid fasta with eUtils.
# 4.  Create customized blast databases.
# 5.  Use magicBlast on one SRRA versus individual plasmid databases.
# 6.  Create individual contig.fasta files and generate a table with the % of plasmid sequences covered by contigs.
# 7.  Generate plasmid and contig visualization using Circos.
# 8.  Create customized contig db and blast contigs reciprocally.
# 9.  Identify contigs matching each other and other plasmids.
# 10. Download and parse plasmid GenBank genome files.
# 11.  BLAST contigs against gapless genome databases.
# 12. Identify contigs matching gapless genome sequences.




### STEPS ###
# ==================

# Step 1.
# --------------------------
      # Extract all genomic gapless genome fasta and gff.

      # from NCBI website, go to Assembly database:
      # staphylococcus aureus[Organism] 
      # Filters: Status = Latest
      # Filters: Assembly level = complete genomes
      # download Assembly IDs
      # extract first column and remove column names
      # save as "gapless_genome_assemblies.txt" in /got_plasmid_gitHub/gapless_genomes/assembly_ID/

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/gapless_genomes/assembly_ID/
      ## chmod +x Path_to/got_plasmid_gitHub/gapless_genomes/assembly_ID/genome_getter.sh
      ## ./genome_getter.sh
      #
      # genome_getter.sh:
      # -----------------
      # #!/bin/bash
      # for i in `cat gapless_genome_assemblies.txt`; do wget `esearch -db assembly -query "$i" | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'`; done
      # gunzip *.gz
      # cat *.fna > gapless_genome_plasmid.fasta
      # mv gapless_genome_plasmid.fasta ../fasta/.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/gapless_genomes/assembly_ID/
      ## chmod +x Path_to/got_plasmid_gitHub/gapless_genomes/assembly_ID/gff_getter.sh
      ## ./gff_getter.sh
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
      ## cd Path_to/got_plasmid_gitHub/

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
      ## cd Path_to/got_plasmid_gitHub/

      ### ./R
      ### use esearch_plasmid.R
      ### q()

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/plasmids/fasta/
      ## chmod +x Path_to/got_plasmid_gitHub/plasmids/fasta/esearch_plasmid.sh
      ## ./esearch_plasmid.sh

# Step 4. 
# ------------------
      # Create customized blast databases.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/gapless_genomes/fasta/
      ## chmod +x PAth_to/got_plasmid_gitHub/plasmids/fasta/makeblastdb_gapless_genomes.sh 
      ## ./makeblastdb_gapless_genomes.sh
      ## cd Path_to/got_plasmid_gitHub/

      ### ./R
      ### use makeblastdb_individual_plasmid.R
      ### q()

      ## ***command line***
      ## cd PAth_to/got_plasmid_gitHub/plasmids/fasta/
      ## chmod +x PAth_to/got_plasmid_gitHub/plasmids/fasta/makeblastdb_individual_plasmid.sh
      ## ./makeblastdb_individual_plasmid.sh

      
# Step 5. 
# ------------------
      # Use magicBlast on one SRRA versus individual plasmid databases.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use magicblast_plasmid.R
      ### q()

      ## ***command line***
      ## cd PAth_to/got_plasmid_gitHub/plasmids/magic_output/
      ## chmod +x PAth_to/got_plasmid_gitHub/plasmids/magic_output/magicblast_plasmid.sh
      ## ./magicblast_plasmid.sh


# Step 6. 
# ------------------
      # Create individual contig.fasta files and generate a table with the % of plasmid sequences covered by contigs.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use Rscript2_magicBlast_output.R
      ### q()

# Step 7. 
# ------------------
      # Generate plasmid and contig visualization using Circos.
      
      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use Rscript3_circos_plasmid_visualization
      ### q()


      ## ***command line***
      ## cd PAth_to/got_plasmid_gitHub/circos/circos-0.69-6/staph_plasmids/conf/
      ## chmod +x Path_to/got_plasmid_gitHub/plasmids/contigs/*.sh
      ## ./perl_local_circos_temp.sh
      ## mv *.png ../../../../outputs/.


# Step 8. 
# ------------------
      # Create customized contig db and blast contigs reciprocally.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use contigs.R
      ### q()

      ## ***command line***
      ## cd PAth_to/got_plasmid_gitHub/plasmids/contigs/
      ## chmod +x Path_to/got_plasmid_gitHub/plasmids/contigs/*.sh
      ## ./Formatblast_contigs.sh
      ## ./reciprocal_contig_Blast.sh
      
      
# Step 9. 
# ------------------
      # Identify contigs matching each other and other plasmids.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use Rscript4_reciprocal_blast_contig.R
      ### q()

      
# Step 10. 
# ------------------
      # Download and parse plasmid GenBank genome files.

      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use Rscript5_GenBank_parser.R
      ### q()
      


# Step 11. 
# ------------------
      # BLAST contigs against gapless genome databases.
      
      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use blast_contigs_vs_gapless_genomes.R
      ### q()
      
      ## ***command line***
      ## cd PAth_to/got_plasmid_gitHub/plasmids/contigs/
      ## chmod +x Path_to/got_plasmid_gitHub/plasmids/contigs/*.sh
      ## ./blast_contigs_vs_gapless_genomes.sh



# Step 12. 
# ------------------
      # Identify contigs matching gapless genome sequences.
    
      ## ***command line***
      ## cd Path_to/got_plasmid_gitHub/
      
      ### ./R
      ### use Rscript6_contig_match_in_genomes.R
      ### q()


# ==============




### DeBugging/Notes
# =================
# if error message:  /bin/bash^M: bad interpreter: No such file or directory 
# vi file.sh type ':set fileformat=unix' wihout quote
# =================
