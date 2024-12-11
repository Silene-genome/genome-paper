#!/bin/bash -l

#SBATCH -A naiss2023-5-113
#SBATCH -p core -n 8
#SBATCH -t 2-00:00:00 
#SBATCH -J variant-calling-haploid --mail-type=FAIL


module load bioinfo-tools
module load freebayes/1.3.2


folder=/crex/proj/uppstore2017241/ddRAD_Silene/

i=X
freebayes -f ${folder}/3-bwa/mapping-Silene-v.4/S.latifolia_v4.0_autosome_x_y_refA.fasta -L ${folder}/4-freebayes/planA_20230831/male-bam-list.txt --gvcf -r chr${i} --report-monomorphic -m 20 -q 20 -E 3 --min-repeat-entropy 1 -V -n 10 --populations ${folder}/4-freebayes/planA_20230831/male-poplist.txt --ploidy=1 > ${folder}/4-freebayes/planA_20230831/chr${i}_haploid-raw-male.gvcf

i=Y
freebayes -f ${folder}/3-bwa/mapping-Silene-v.4/S.latifolia_v4.0_autosome_x_y_refA.fasta -L ${folder}/4-freebayes/planA_20230831/male-bam-list.txt --gvcf -r chr${i} --report-monomorphic -m 20 -q 20 -E 3 --min-repeat-entropy 1 -V -n 10 --populations ${folder}/4-freebayes/planA_20230831/male-poplist.txt --ploidy=1 > ${folder}/4-freebayes/planA_20230831/chr${i}_haploid-raw-male.gvcf
