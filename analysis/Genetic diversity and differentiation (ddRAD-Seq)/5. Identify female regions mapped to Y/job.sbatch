#!/bin/bash -l

#SBATCH -A naiss2023-5-113  
#SBATCH -p core  -n 2                                                                                                                                          
#SBATCH -t 20:00:00  

module load bioinfo-tools
module load samtools
module load BEDTools/2.29.2

outdir=/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/result_new_20230904/region_female_mapped_Y

samtools depth -a -f ${outdir}/female_goodDP_planA.txt -b ${outdir}/chrY.region.bed -H -o ${outdir}/chrY.DP.tsv

perl filter_dp.pl chrY.DP.tsv 3 > chrY.minDP3.count

awk '$3>=2 {print $0}' chrY.minDP3.count > chrY.minDP3.minIND2.txt

awk '{print $1"\t"$2-1"\t"$2}' chrY.minDP3.minIND2.txt > chrY.minDP3.minIND2.raw.bed

sort -k1,1 -k2,2n chrY.minDP3.minIND2.raw.bed > chrY.minDP3.minIND2.sort.bed

bedtools merge -i ./chrY.minDP3.minIND2.sort.bed  > chrY.minDP3.minIND2.merge.bed
