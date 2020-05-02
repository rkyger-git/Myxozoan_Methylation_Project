## C-shasta-WGBS-qc-and-bss-analysis.md

#### This document shows the linux commands I used for sequence quality control, alignment, and analysis of bisulfite sequence data for C. shasta

### Sequence Quality Control with FastQC
```
#Input files:
1B-C-shasta_R1.fastq.gz #forward reads
1B-C-shasta_R2.fastq.gz #reverse reads

#run fastqc on forward and reverse WGBS reads for Cshasta:
fastqc 1B-C-shasta_R1.fastq.gz
fastqc 1B-C-shasta_R2.fastq.gz

#run trim_galore on forward and reverse WGBS reads, “--paired” specifies that the reads are paired end, “--trim1” trims one additional base pair from the 3' end of both reads, 
#this is necessary because Bowtie 1 considers self-containing sequences to be invalid paired-end alignments:
trim_galore --paired --trim1 1B-C-shasta_R1.fastq.gz 1B-C-shasta_R2.fastq.gz

#Output files:
1B-C-shasta_R1_val_1.fq.gz #forward reads
1B-C-shasta_R2_val_2.fq.gz #reverse reads

#run fastqc on trimmed forward and reverse WGBS reads for Cshasta:
fastqc 1B-C-shasta_R1_val_1.fq.gz
fastqc 1B-C-shasta_R2_val_2.fq.gz
```

### WGBS Sequence Alignment using Bismark with Bowtie 1
```
#Note: for each of the below: “--bowtie1” tells bismark to use bowtie 1.

#create bisulfite reference genome for lambda phage, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . enterobacteria-phage-lambda-genome.fasta 

#align reads to bisulfite reference genome, “--un” tells bismark to keep all reads that don’t map to the reference genome, “--gzip” tells bimark to compress the unmapped reads, 
#“-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --un --gzip -n 1 . -1 1B-C-shasta_R1_val_1.fq.gz -2 1B-C-shasta_R2_val_2.fq.gz

#BAM output file:
1B-C-shasta_R1_val_1_bismark_pe.bam #phage alignment


#create bisulfite reference genome for rainbow trout, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . GCA_002163495.1_Omyk_1.0_genomic.fa

#align unmapped reads to bisulfite reference genome, “--un” tells bismark to keep all reads that don’t map to the reference genome, “--gzip” tells bimark to compress the unmapped reads, 
#“-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --un --gzip -n 1 . -1 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz -2 1B-C-shasta_R2_val_2.fq.gz_unmapped_reads_2.fq.gz

#BAM output file:
1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.bam #trout alignment 


#create bisulfite reference genome for C.shasta, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . cleaned_genome_Cshasta.fasta

#align unmapped reads to bisulfite reference genome,“--gzip” tells bimark to compress the unmapped reads, “-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),
#“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --gzip -n 1 . -1 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz -2 1B-C-shasta_R2_val_2.fq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz

#BAM output file:
1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.bam #Cshasta alignment
```

### Remove PCR duplication bias
```
#remove PCR duplication bias from the BAM output files for produced in the previous three steps, “--bam” tells bismark to output a BAM file (instead of SAM): 
perl deduplicate_bismark --bam 1B-C-shasta_R1_val_1_bismark_pe.bam #phage alignment
perl deduplicate_bismark --bam 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.bam #trout alignment 
perl deduplicate_bismark --bam 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.bam #Cshasta alignment

#Deduplicated BAM output files:
1B-C-shasta_R1_val_1_bismark_pe.deduplicated.sam #phage alignment
1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #trout alignment
1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #Cshasta alignment
```

### Acquire the deduplicated methylation data from deduplicated BAM files
```
#acquire the deduplicated methylation data from deduplicated BAM files, “-p” tells bismark that the data is paired-end, 
#“--no_overlap” tells bismark to avoid double counting of nucleotides for overlapping reads, “--gzip” tells bimark to compress all output files:
perl bismark_methylation_extractor -p --no_overlap --gzip 1B-C-shasta_R1_val_1_bismark_pe.deduplicated.bam #phage alignment
perl bismark_methylation_extractor -p --no_overlap --gzip 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #trout alignment
perl bismark_methylation_extractor -p --no_overlap --gzip 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #Cshasta alignment
```

### Acquire the nucleotide coverage
```
#acquire the nucleotide coverage from deduplicated BAM file for the final C.shasta alignment, used “--genome_folder” to specify the full path to bisulfite reference genome for C.shasta:
perl bam2nuc --genome_folder /srv/WolfCloudNFS/HomeWolfCloud1046/rkyger/myxozoa/cshasta-bs/cshasta-bs-clean/ 1B-C-shasta_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam
```

### Create HTML reports
```
#This takes report files produced in previous steps and creates an html report file. 
#I ran this command three times (phage data, trout data, C.shasta data)
perl bismark2report

#renamed html report files:
WGBS-C-shasta-aligned-to-phage-bismark-report.html #phage data
WGBS-C-shasta-aligned-trout-bismark-report.html #trout data
WGBS-C-shasta-aligned-to-C-shasta-bismark-report.html #C.shasta data
```
