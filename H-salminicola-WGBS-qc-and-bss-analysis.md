## H-salminicola-WGBS-qc-and-bss-analysis.md

#### This document shows the linux commands used for sequence quality control, alignment, and analysis of bisulfite sequence data for H. saliminicola

### Sequence Quality Control with FastQC
```
#Input files:
2B-H-salmonicola_R1.fastq.gz #forward reads
2B-H-salmonicola_R2.fastq.gz #reverse reads

#run fastqc on forward and reverse WGBS reads for H. saliminicola:
fastqc 2B-H-salmonicola_R1.fastq.gz
fastqc 2B-H-salmonicola_R2.fastq.gz

#run trim_galore on forward and reverse WGBS reads, “--paired” specifies that the reads are paired end, “--trim1” trims one additional base pair from the 3' end of both reads, 
#this is necessary because Bowtie 1 considers self-containing sequences to be invalid paired-end alignments:
trim_galore --paired --trim1 2B-H-salmonicola_R1.fastq.gz 2B-H-salmonicola_R2.fastq.gz

#Output files:
2B-H-salmonicola_R1_val_1.fq.gz #forward reads
2B-H-salmonicola_R1_val_2.fq.gz #reverse reads

#run fastqc on trimmed forward and reverse WGBS reads for Hsaliminicola:
fastqc 2B-H-salmonicola_R1_val_1.fq.gz
fastqc 2B-H-salmonicola_R1_val_2.fq.gz
```

### WGBS Sequence Alignment using Bismark with Bowtie 1
```
#Note: for each of the below: “--bowtie1” tells bismark to use bowtie 1.

#create bisulfite reference genome for lambda phage, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . enterobacteria-phage-lambda-genome.fasta 

#align reads to bisulfite reference genome, “--un” tells bismark to keep all reads that don’t map to the reference genome, “--gzip” tells bimark to compress the unmapped reads, 
#“-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --un --gzip -n 1 . -1 2B-H-salmonicola_R1_val_1.fq.gz -2 2B-H-salmonicola_R1_val_2.fq.gz

#BAM output file:
2B-H-salmonicola_R1_val_1_bismark_pe.bam #phage alignment


#create bisulfite reference genome for chinook salmon, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . GCA_002872995.1_Otsh_v1.0_genomic.fa

#align unmapped reads to bisulfite reference genome, “--un” tells bismark to keep all reads that don’t map to the reference genome, “--gzip” tells bimark to compress the unmapped reads, 
#“-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --un --gzip -n 1 . -1 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz -2 2B-H-salmonicola_R2_val_2.fq.gz_unmapped_reads_2.fq.gz

#BAM output file:
2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.bam #salmon alignment


#create bisulfite reference genome for H.saliminicola, “.” tells bismark to put the new bisulfite reference genome in the current directory:
perl bismark_genome_preparation --bowtie1 . IDBA_HZS_AllfilterFinal2_300.fasta

#align unmapped reads to bisulfite reference genome,“--gzip” tells bimark to compress the unmapped reads, “-n 1” tells bimsark to tolerate one non-bisulfite mismatch per read (recommended setting),
#“.” tells bismark that the reference genome is in the current directory, “-1” and “-2” specify the forward and reverse reads:
perl bismark --bowtie1 --gzip -n 1 . -1 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1.fq.gz -2 2B-H-salmonicola_R2_val_2.fq.gz_unmapped_reads_2.fq.gz_unmapped_reads_2.fq.gz

#BAM output file:
2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.bam #H.salminicola alignment
```

### Remove PCR duplication bias
```
#remove PCR duplication bias from the BAM output files for produced in the previous three steps, “--bam” tells bismark to output a BAM file (instead of SAM): 
perl deduplicate_bismark --bam 2B-H-salmonicola_R1_val_1_bismark_pe.bam #phage alignment
perl deduplicate_bismark --bam 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.bam #salmon alignment
perl deduplicate_bismark --bam 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.bam #H.salminicola alignment

#Deduplicated BAM output files:
2B-H-salmonicola_R1_val_1_bismark_pe.deduplicated.bam #phage alignment
2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #salmon alignment
2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #H.salminicola alignment
```

### Acquire the deduplicated methylation data from deduplicated BAM files
```
#acquire the deduplicated methylation data from deduplicated BAM files, “-p” tells bismark that the data is paired-end, 
#“--no_overlap” tells bismark to avoid double counting of nucleotides for overlapping reads, “--gzip” tells bimark to compress all output files:
perl bismark_methylation_extractor -p --no_overlap 2B-H-salmonicola_R1_val_1_bismark_pe.deduplicated.bam #phage alignment
perl bismark_methylation_extractor -p --no_overlap 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #salmon alignment
perl bismark_methylation_extractor -p --no_overlap --gzip 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam #H.salminicola alignment
```

### Acquire the nucleotide coverage
```
#acquire the nucleotide coverage from deduplicated BAM file for the final H.salminicola alignment, used “--genome_folder” to specify the full path to bisulfite reference genome for H.saliminicola
perl bam2nuc --genome_folder /srv/WolfCloudNFS/HomeWolfCloud1046/rkyger/myxozoa/hsalminicola-bs/hsaliminicola-bs-clean 2B-H-salmonicola_R1_val_1.fq.gz_unmapped_reads_1.fq.gz_unmapped_reads_1_bismark_pe.deduplicated.bam
```

### Create HTML reports
```
#This takes report files produced in previous steps and creates an html report file. 
#I ran this command three times (phage data, trout data, H.salmincola data)
perl bismark2report

#renamed html report files:
WGBS-H-salmonicola-aligned-to-phage-bismark-report #phage data
WGBS-H-salmonicola-aligned-to-salmon-bismark-report #trout data
WGBS-H-salmonicola-aligned-to-H-salminicola-bismark-report #H.salminicola data
```
