#!/bin/bash -l
#SBATCH -A naiss2023-22-242
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J index

cd /crex/proj/uppstore2017241/ddRAD_Silene/3-bwa/mapping-Silene-v.4/

module load bioinfo-tools bwa/0.7.17

#bwa index S.latifolia_v4.0_autosome_x_y_refA.fasta

bwa index S.latifolia_v4.0_autosome_x_refB.fasta

