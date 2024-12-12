

# the Silene latofolia paper

All the codes to generate the assembly, annotation and analysis described in our manuscript "The Silene latifolia genome and its giant Y chromosome".


Abstract
============

In some species, the Y is a tiny chromosome but the dioecious plant Silene latifolia has a giant ∼550 Mb Y chromosome, which has remained unsequenced so far. Here we used a hybrid approach to obtain a high-quality male S. latifolia genome. Using mutants for sexual phenotype, we identified candidate sex-determining genes on the Y. Comparative analysis of the sex chromosomes with outgroups showed the Y is surprisingly rearranged and degenerated for a ∼11 MY-old system. Recombination suppression between X and Y extended in a stepwise process, and triggered a massive accumulation of repeats on the Y, as well as in the non-recombining pericentromeric region of the X, leading to giant sex chromosomes.


## Assembly

1. First, we run Flye 2.8.1-b1676:

```
Flye -t 64 --nano-raw ${LONG_READS} -g 3g -o flye-silene-ONT-100X --asm-coverage 100 --min-overlap 10000

```

2. Polishing, steps 1 and 2 using minimap version 0.1.0, RACON v1.4.9, hthit version 0.1.0 and ntedit version 1.2.3:

```
.DELETE_ON_ERROR:
#racon version v1.4.9
RACON=/miniconda/envs/polishing/bin/racon
#minimap2 version 2.15-r905
MM=/miniconda/envs/polishing/bin/minimap2
# nthits version 0.1.0
NTH=/miniconda/envs/polishing/bin/nthits
# ntedit version 1.2.3
NTE=/miniconda/envs/polishing/bin/ntedit
#TIME=/path/time
#LONG-READ POLISHING
#first round of long-read polishing
$(PREF).r1.paf:
        ${MM} -x map-ont -t ${CPU} ${ASM} ${READS} > $@
$(PREF).r1.fa:$(PREF).r1.paf
        ${RACON} -u -t ${CPU} ${READS} $< ${ASM} > $@ 2> $(PREF).r1.racon.log
#second round of long-read polishing
$(PREF).r2.paf:$(PREF).r1.fa
        ${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r2.fa:$(PREF).r2.paf
        ${RACON} -t ${CPU} ${READS} $< $(PREF).r1.fa > $@ 2> $(PREF).r2.racon.log
#3th round of long-read polishing
$(PREF).r3.paf:$(PREF).r2.fa
        ${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r3.fa:$(PREF).r3.paf
        ${RACON} -t ${CPU} ${READS} $< $(PREF).r2.fa > $@ 2> $(PREF).r3.racon.log
#4th round of long-read polishing
$(PREF).r4.paf:$(PREF).r3.fa
        ${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r4.fa:$(PREF).r4.paf
        ${RACON} -t ${CPU} ${READS} $< $(PREF).r3.fa > $@ 2> $(PREF).r4.racon.log
#5th round of long-read polishing
$(PREF).r5.paf:$(PREF).r4.fa
        ${MM} -x map-ont -t ${CPU} $< ${READS} > $@
$(PREF).r5.fa:$(PREF).r5.paf
        ${RACON} -t ${CPU} ${READS} $< $(PREF).r4.fa > $@ 2> $(PREF).r5.racon.log
#SHORT-READ POLISHING
$(PREF).nthits_k60.bf:$(PREF).r5.fa
        echo ${SREADS} | xargs -n 1 > shortreads.txt
        ${NTH} -b 36 -k 40 -t ${CPU} -p $(PREF).nthits --outbloom --solid @shortreads.txt
        ${NTH} -b 36 -k 50 -t ${CPU} -p $(PREF).nthits --outbloom --solid @shortreads.txt
        ${NTH} -b 36 -k 60 -t ${CPU} -p $(PREF).nthits --outbloom --solid @shortreads.txt
# three rounds of short-read polishing with ntEdit after two round of RACON
$(PREF).r5.nt3_edited.fa:$(PREF).r5.fa $(PREF).nthits_k60.bf
        ${NTE} -t ${CPU} -i 5 -d 5 -b $(PREF).r5.nt1 -r $(PREF).nthits_k60.bf -f $<
        ${NTE} -t ${CPU} -i 5 -d 5 -b $(PREF).r5.nt2 -r $(PREF).nthits_k50.bf -f $(PREF).r5.nt1_edited.fa
        ${NTE} -t ${CPU} -i 5 -d 5 -b $(PREF).r5.nt3 -r $(PREF).nthits_k40.bf -f $(PREF).r5.nt2_edited.fa
# all done
all: $(PREF).r5.nt3_edited.fa


``` 

2. Polishing, step 3 using pilon:

```
#parameters file
cat params-polishing-v01.yml 
fasta : asm/flye-polish-100-v2.0.polished.fa
reads : 'reads/CSR_ABOSDE_1_{1,2}_HGHFNDSXY.12BA301_clean.fastq.gz'
nsplit : 2000000

cat params-polishing-v02.yml 
fasta : pilon_round1/flye-polish-100-v2.0.pilon_round1.fa
reads : 'reads/CSR_ABOSDE_1_{1,2}_HGHFNDSXY.12BA301_clean.fastq.gz'
nsplit : 2000000
cpu : 40
outdir : pilon_round2
mem_mem : 80

cat params-polishing-v03.yml 
fasta : pilon_round2/flye-polish-100-v2.0.pilon_round2.fa
reads : 'reads/CSR_ABOSDE_1_{1,2}_HGHFNDSXY.12BA301_clean.fastq.gz'
nsplit : 2000000
cpu : 40
outdir : pilon_round3
mem_mem : 80

#nextflow commands
#round 1
nextflow run ppilon/main.nf -profile singularity -c nextflow.leftraru.config -params-file params-polishing-v01.yml 

#round 2
nextflow run ppilon/main.nf -profile singularity -c nextflow.leftraru.config -params-file params-polishing-v02.yml 

#round 3
nextflow run ppilon/main.nf -profile singularity -c nextflow.leftraru.config -params-file params-polishing-v03.yml 


```

3. To compute metrics we run samtools and a script in perl:

```
samtools faidx $1.fa
perl compute-stats-from-fai.pl -a $1.fa.fai

```

4. we run BUSCO analysis using singularity:

```
export PATH=/beegfs/data/soft/singularity-3.1.1/bin/:$PATH
CONTAINER=busco_v5.0.0_cv1.sif

singularity exec -B /beegfs:/beegfs  busco -i fasta-file -l viridiplantae_odb10 -o anotation-eugene.BUSCO -c 20 -m prot

```

5. Genome size reduction, we run windowmasker and purge_dups:

```
windowmasker -mk_counts -in flye-polish-100-v2.0.pilon_round3.fa -out flye-polish-100-v2.0.pilon_round3.windowmask.out
windowmasker -ustat flye-polish-100-v2.0.pilon_round3.windowmask.out -in flye-polish-100-v2.0.pilon_round3.fa -out flye-polish-100-v2.0.pilon_round3.windowmask.step2.out -dust true

```

```
.DELETE_ON_ERROR:


ASM=flye-polish-100-v2.0.pilon_round3.fa
R1=ont_silene_raw.fastq.gz
R2=ont_silene_raw.r2.fastq.gz
PURGE=purge_haplotigs
#map the reads long reads
CPU=40
ont_silene.bam:
        -rm -f long_reads.fq.gz
        mkfifo long_reads.fq.gz
        cat ${R2} ${R1} > long_reads.fq.gz &
        minimap2 -t ${CPU} -ax map-ont ${ASM} long_reads.fq.gz --secondary=no | samtools view -b - | samtools sort -@10 -m 2G -o $@ -T tmp.aln
        samtools index $@
#compute coverage with purge_dups
ont_silene.bam.gencov:ont_silene.bam
        purge_haplotigs  hist  -b $<  -g ${ASM} -t ${CPU} > silene_coverage_ont.txt
#loking at the plot we do define l=5, m=30, h=60
silene_coverage_stats.csv:ont_silene.bam.gencov
        purge_haplotigs cov  -i $< -l 5 -m 30 -h 60  -o $@
#we purge the haplotypes without using the repeats
silene_purge_r1.fasta:silene_coverage_stats.csv
        ${PURGE}  purge  -t ${CPU} -g ${ASM}  -c $< -o silene_purge_r1 -d -b ont_silene.bam > silene_purge_r1.log
#we purge the haplotypes using the window masker
silene_purge_r2.fasta:silene_coverage_stats.csv
        ${PURGE}  purge  -t ${CPU} -g ${ASM} -r flye.windowmask.bed -c $< -o silene_purge_r2 -d -b ont_silene.bam > silene_purge_r2.log
#we seq all the passes
all: ont_silene.bam ont_silene.bam.gencov silene_coverage_stats.csv silene_purge_r1.fasta silene_purge_r2.fasta

```


### BrumiR outputs

| File  |  Description  |   
|-------|---------------|
| prefix.brumir.candidate_miRNA.fasta   |  fasta file with all the candidates with their KM and KC values respectively. |
|  prefix.brumir.other_sequences.txt |  asta file with all long sequences expressed in the sample, they are putative long non-coding RNAs. |
| prefix.brumir.RFAM_HITS.txt | table with a list of putative tRNAs or rRNAs present in the RFAM database. |

 