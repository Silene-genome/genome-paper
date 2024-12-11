### RNASeq data analysis ###

# reference genome were indexed using STAR

STAR --runMode genomeGenerate  --genomeFastaFiles  genome.fasta --sjdbGTFfile genome.gtf --genomeDir genome_index/ 

# Samples were process separately using the following command lines, here as an exemple sample M-S5-R1 (male stade 5, replicate 1):
# As described in the M&M, female and male were aligned on different reference (respectively without and with chrY)
cutadapt  -a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC  -G AGATCGGAAGAGCACACGT -q 20 --minimum-length 50 \ 
 -o trimmed/M-S5-R1.1.fastq -p trimmed/M-S5-R1.2.fastq data/M-S5-R1.1.fastq  data/M-S5-R1.2.fastq

STAR --genomeDir genome_index/ --readFilesIn trimmed/M-S5-R1.1.fastq trimmed/M-S5-R1.2.fastq --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterType BySJout --outTmpDir STARtmp --outFileNamePrefix M-S5-R1/ 

samtools sort -n -o M-S5-R1/M-S5.sorted1.bam  M-S5-R1/Aligned.sortedByCoord.out.bam

samtools fixmate -m M-S5-R1/M-S5-R1.sorted1.bam M-S5-R1/M-S5-R1.sorted1.fixmate.bam

samtools sort M-S5-R1/M-S5-R1.sorted1.fixmate.bam M-S5-R1/M-S5-R1.sorted2.fixmate.bam

samtools view -F 1548 -b -q 50 M-S5-R1/M-S5-R1.sorted2.fixmate.bam > M-S5-R1/M-S5-R1.sorted.nodup.filt.bam

samtools index M-S5-R1/M-S5-R1.sorted.nodup.filt.bam

# As detailled in the M&M FeatureCcounts was used to obtain the count table

featureCounts -s 0 -t exon  -g gene_id -p -F GTF -a genome.gtf -G genome.fasta -o featureCounts/count_table M-S5-R1/M-S5-R1.sorted.nodup.filt.bam M-S5-R2/M-S5-R2.sorted.nodup.filt.bam M-S5-R3/M-S5-R3.sorted.nodup.filt.bam M-S8-R1/M-S8-R1.sorted.nodup.filt.bam M-S8-R2/M-S8-R2.sorted.nodup.filt.bam M-S8-R3/M-S8-R3.sorted.nodup.filt.bam F-S5-R1/F-S5-R1.sorted.nodup.filt.bam F-S5-R2/F-S5-R2.sorted.nodup.filt.bam F-S5-R3/F-S5-R3.sorted.nodup.filt.bam  F-S8-R1/F-S8-R1.sorted.nodup.filt.bam F-S8-R2/F-S8-R2.sorted.nodup.filt.bam F-S8-R3/F-S8-R3.sorted.nodup.filt.bam

# Reformat Feature count output to be able to use it for differential expression analysis: Get ride of first line (command line), and rename columns. 
# Create Metadata file. First column = name of sample identical to the one in the count table. The columns after contain the metadata.
# everything put in DE_analysis/input

# Run Differential expression analysis using R
cd DE_analysis
cd Quality_Control
# Quality control script need a working installation of https://forgemia.inra.fr/GNet/dicoexpress/
Rscript Script_quality_control.R
# Run deseq2
cd ..
mkdir deseq2_analysis
Rscript Script_deseq2.R