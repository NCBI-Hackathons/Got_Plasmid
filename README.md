# Got Plasmids?
## (or phages, repetitive elements or polintons?)

## Goal
Defining the clinically relevant plasmids of Staphylococcus aureus and other pathogenic bacteria.
Challenge: mobile genetic elements, including plasmids and phages,  are ambiguous as rich in repeats and low complexity sequences.

AIMS:
- Identification of unique plasmids in individual patient samples
- Differentiation between screening (nose swabs) and lesions
- Building a rapid pipeline to differentiate between plasmids in these differential sites.  
- Extending that to other bacterial species

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

## Dependencies:
* R
* R packages (data.table, rowr, GenomicRanges, IRanges, stringr, ggplot2)
* Blast
* magicBlast
* eUtils


## WorkFlow
![My image](https://github.com/NCBI-Hackathons/Pathogenic_Pangenomes/blob/master/images/workflow_2.png)
