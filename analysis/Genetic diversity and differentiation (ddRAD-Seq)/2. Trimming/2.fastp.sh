mkdir ./ddRAD_Silene/2-fastp/

#copy path to the fasta of the selected populations: ESP1, POL1, POL3
find ./ddRAD_Silene/1-demultiplex/Sample_RW*/ESP1* > ./ddRAD_Silene/2-fastp/samples_list.txt
find ./ddRAD_Silene/1-demultiplex/Sample_RW*/POL1* >> ./ddRAD_Silene/2-fastp/samples_list.txt
find ./ddRAD_Silene/1-demultiplex/Sample_RW*/POL3* >> ./ddRAD_Silene/2-fastp/samples_list.txt

find ./ddRAD_Silene/1-demultiplex/Sample_RW*/BEL1* > ./ddRAD_Silene/2-fastp/samples_list2.txt
find ./ddRAD_Silene/1-demultiplex/Sample_RW*/FRA1* >> ./ddRAD_Silene/2-fastp/samples_list2.txt

#remove rem files from the list
sed '/rem/d' ./ddRAD_Silene/2-fastp/samples_list.txt > ./ddRAD_Silene/2-fastp/fasta_list.txt
#separate forward - reverse from the list
sed '/1\./!d' ./ddRAD_Silene/2-fastp/fasta_list.txt > ./ddRAD_Silene/2-fastp/fasta-r1_list.txt
sed '/2\./!d' ./ddRAD_Silene/2-fastp/fasta_list.txt > ./ddRAD_Silene/2-fastp/fasta-r2_list.txt
# concatenate in column r1 and r2 with one sample per line, use semi-column to delimitate the column
paste ./ddRAD_Silene/2-fastp/fasta-r1_list.txt ./ddRAD_Silene/2-fastp/fasta-r2_list.txt > ./ddRAD_Silene/2-fastp/fasta-r.txt
#add a column with sample names
#extract sample names from fasta-r.txt file (did it with notepad)
paste ./ddRAD_Silene/2-fastp/fasta-r.txt ./ddRAD_Silene/2-fastp/sample-names.txt > ./ddRAD_Silene/2-fastp/fasta-r-names.txt

#trimming with fastp
mkdir ./ddRAD_Silene/0-scripts/fastp/

list=fasta-r-names

IFS=';'
while read line
do
linearray=( $line )
r1=${linearray[0]}
r2=${linearray[1]}
samplename=${linearray[2]}

echo '#!/bin/bash -l' > ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo '#SBATCH -A naiss2023-22-242' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo '#SBATCH -p core' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo '#SBATCH -n 2' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo '#SBATCH -t 01:00:00' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo '#SBATCH -J fastp --mail-type=FAIL,TIME_LIMIT_80' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo 'cd /crex/proj/uppstore2017241/' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo 'module load bioinfo-tools fastp/0.23.2' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo 'outdir=./ddRAD_Silene/2-fastp' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo 'TW="--length_required 60" # half of read length' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh
echo 'fastp -i '$r1' -I '$r2' -o $outdir/'$samplename'.R1.fq.gz -O $outdir/'$samplename'.R2.fq.gz --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction $TW -q 15 -u 50 -j $outdir/'$samplename'.json -h $outdir/'$samplename'.html --detect_adapter_for_pe &> $outdir/'$samplename'.trim.log ' >> ./ddRAD_Silene/0-scripts/fastp/$samplename.sh

sbatch ./ddRAD_Silene/0-scripts/fastp/$samplename.sh

done < ./ddRAD_Silene/2-fastp/fasta-r-names.txt

#per read cutting by quality options
  # -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
  # -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
  # -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
  # -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
  # -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
#base correction by overlap analysis options
#-c, --correction                   enable base correction in overlapped regions (only for PE data), default is disabled
  # --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
  # --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
  # --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
#-q, --qualified_quality_phred      the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
#-u unqualified percent limit, default 40
# reporting options
#-j, --json                         the json format report file name (string [=fastp.json])
#-h, --html                         the html format report file name (string [=fastp.html])
#detect_adapter_for_pe          by default, the adapter sequence auto-detection is enabled for SE data only, turn on this option to enable it for PE data.
    
#-y Low complexity filter The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
#its default value is 30
