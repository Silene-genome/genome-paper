##############################################################
# Loop for individual alignment on a ref.
##############################################################

#create directory
mkdir ./ddRAD_Silene/3-bwa/
mkdir ./ddRAD_Silene/0-scripts/bwa/

#decompose the genome - see bash script
#create the new reference
cd ./mapping-Silene-v.4/decomposed-genome/
cat chr1.fasta chr2.fasta chr3.fasta chr4.fasta chr5.fasta chr6.fasta chr7.fasta chr8.fasta chr9.fasta chr10.fasta chr11.fasta chrX.fasta chrY.fasta > S.latifolia_v.4.0_autosome_x_y_refA.fasta
#minus y for refB

#software to load
module load bioinfo-tools bwa/0.7.17

#create the INDEX - use the bash script
#build index
bwa index REF

#list
find ./ddRAD_Silene/2-fastp/*.fq.gz > ./ddRAD_Silene/3-bwa/list2.txt
#separate forward - reverse from the list
sed '/R1\./!d' ./ddRAD_Silene/3-bwa/list2.txt > ./ddRAD_Silene/3-bwa/listR1.txt
sed '/R2\./!d' ./ddRAD_Silene/3-bwa/list2.txt > ./ddRAD_Silene/3-bwa/listR2.txt
# concatenate in column r1 and r2 with one sample per line, use semi-column to delimitate the column
paste -d ';' ./ddRAD_Silene/3-bwa/listR1.txt ./ddRAD_Silene/3-bwa/listR2.txt > ./ddRAD_Silene/3-bwa/listR.txt
#add a column with sample names
#extract sample names from fasta-r.txt file (did it with notepad)
paste -d ';' ./ddRAD_Silene/3-bwa/listR.txt ./ddRAD_Silene/3-bwa/names.txt > ./ddRAD_Silene/3-bwa/listR-names.txt>

list=female-list.txt
list=male-list.txt

#output directories
mkdir ./planA/
mkdir ./planB/


IFS=';'
while read line
do 
linearray=( $line )
r1=${linearray[0]}
r2=${linearray[1]}
sample=${linearray[2]}

echo '#!/bin/bash -l' > ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo '#SBATCH -A naiss2023-22-242' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo '#SBATCH -p core' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo '#SBATCH -n 8' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo '#SBATCH -t 10:00:00' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo '#SBATCH -J bwamem --mail-type=FAIL' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'cd /crex/proj/uppstore2017241/' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'REF=S.latifolia_v4.0_autosome_x_y_refA.fasta' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'indir=./ddRAD_Silene/3-bwa/mapping-Silene-v.4' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'outdir=./ddRAD_Silene/3-bwa/mapping-Silene-v.4/planA' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'individual='$sample'' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'n=8' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'module load bioinfo-tools bwa/0.7.17' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'module load bioinfo-tools samtools/1.16' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'bwa mem -t $n -R "@RG\tID:"$individual"p\tSM:"$individual"\tPL:ILLUMINA\tPI:330" $indir/$REF '$r1' '$r2' | samtools view -bS | samtools sort -o $outdir/$individual.bam' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'samtools index $outdir/$individual.bam' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh
echo 'samtools stat $outdir/$individual.bam > $outdir/$individual-stats.txt' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh 
echo 'samtools coverage $outdir/$individual.bam > $outdir/$individual-coverage.txt' >> ./ddRAD_Silene/0-scripts/bwa/$sample.sh

sbatch ./ddRAD_Silene/0-scripts/bwa/$sample.sh

done < ./ddRAD_Silene/3-bwa/mapping-Silene-v.4/female-list.txt

#-R STR: complete the read group header line. It is required by GATK. but did not work with vimnalis data
#index bam after marking duplicates with picard

#For chromosome Y planB see bwa-Y.sh script
