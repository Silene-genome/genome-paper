#!/bin/bash -l
#SBATCH -A naiss2023-22-242
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J variant-calling --mail-type=FAIL,TIME_LIMIT_80

#run per chrm with all the samples at once
#exclude the samples with little data

#Loop for chromosome 1 to 11, and X for females

# males
# for i in {1..11}
# do
    # echo $i
    # echo chr${i}
	# (echo '#!/bin/bash -l';
	 # echo "module load bioinfo-tools";
	 # echo "module load freebayes/1.3.2";
	 # echo "freebayes -f ./ddRAD_Silene/3-bwa/mapping-Silene-v.4/S.latifolia_v4.0_autosome_x_y_refA.fasta -L ./ddRAD_Silene/4-freebayes/planA_20230831/male-bam-list.txt --gvcf -r chr${i} --report-monomorphic -m 20 -q 20 -E 3 --min-repeat-entropy 1 -V -n 10 --populations ./ddRAD_Silene/4-freebayes/planA_20230831/male-poplist.txt > ./ddRAD_Silene/4-freebayes/planA_20230831/chr${i}-raw-male.gvcf" ) | sbatch -A naiss2023-22-242 -p core -n 10 -t 2-00:00:00 -J variant-call --mail-type=FAIL,TIME_LIMIT_80

# done

#turn off/ remove "--genotype-qualities"

# females

# for i in {1..11}
# do
    # echo $i
    # echo chr${i}
	# (echo '#!/bin/bash -l';
	 # echo "module load bioinfo-tools";
	 # echo "module load freebayes/1.3.2";
	 # echo "freebayes -f ./ddRAD_Silene/3-bwa/mapping-Silene-v.4/S.latifolia_v4.0_autosome_x_refB.fasta -L ./ddRAD_Silene/4-freebayes/planB_20230831/female-bam-list.txt --gvcf -r chr${i} --report-monomorphic -m 20 -q 20 -E 3 --min-repeat-entropy 1 -V -n 10 --populations ./ddRAD_Silene/4-freebayes/planB_20230831/female-poplist.txt > ./ddRAD_Silene/4-freebayes/planB_20230831/chr${i}-raw-female.gvcf" ) | sbatch -A naiss2023-22-242 -p core -n 10 -t 2-00:00:00 -J variant-call --mail-type=FAIL,TIME_LIMIT_80

# done

for i in X
do
    echo $i
    echo chr${i}
	(echo '#!/bin/bash -l';
	 echo "module load bioinfo-tools";
	 echo "module load freebayes/1.3.2";
	 
	 echo "freebayes -f ./ddRAD_Silene/3-bwa/mapping-Silene-v.4/S.latifolia_v4.0_autosome_x_refB.fasta -L ./ddRAD_Silene/4-freebayes/planB_20230831/female-bam-list.txt --gvcf -r chr${i} --report-monomorphic -m 20 -q 20 -E 3 --min-repeat-entropy 1 -V -n 10 --populations ./ddRAD_Silene/4-freebayes/planB_20230831/female-poplist.txt > ./ddRAD_Silene/4-freebayes/planB_20230831/chr${i}-raw-female.gvcf" ) | sbatch -A naiss2023-22-242 -p core -n 10 -t 2-00:00:00 -J variant-call --mail-type=FAIL,TIME_LIMIT_80

done

#All at once

# for i in $(seq 1 11) X
# do 
# echo $i
# done
