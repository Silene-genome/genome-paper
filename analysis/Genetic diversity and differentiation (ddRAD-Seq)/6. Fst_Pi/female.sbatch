#!/bin/bash -l

#SBATCH -p devcore -n 1                                                                                                                                          
#SBATCH -t 1:00:00                                                                                                                                            
#SBATCH -A naiss2023-22-242                                                                                                                                    

module load bioinfo-tools
module load bcftools
module load vcftools
module load vcflib
module load pixy


chr=$1
indir=/proj/uppstore2017241/ddRAD_Silene/4-freebayes/planB_120623/female_120623
input=${indir}/chr${chr}-raw.gvcf
outdir=/crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/result_new/female
output=${outdir}/chr${chr}_norm.gvcf


# step 1 norm
vcfallelicprimitives -kg ${input} > ${output} 

# step 2 filter, minDP of 6, remove indel, keep target individuals
vcftools --vcf ${outdir}/chr${chr}_norm.gvcf \
  --minDP 6 \
  --remove-indels \
  --recode-INFO-all \
  --keep ${outdir}/sample_goodDP.female.tsv \
  --recode --stdout | bgzip -c > ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz

# # step 3 filter data missingness per pop
 vcftools --gzvcf  ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz    --keep ${outdir}/BEL1_goodDP.list --missing-site --out ${outdir}/BEL1.chr${chr}
 vcftools --gzvcf  ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz    --keep ${outdir}/FRA1_goodDP.list --missing-site --out ${outdir}/FRA1.chr${chr}
 vcftools --gzvcf  ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz    --keep ${outdir}/POL3_goodDP.list --missing-site --out ${outdir}/POL3.chr${chr}


 thres=0.5
 cat ${outdir}/BEL1.chr${chr}.lmiss ${outdir}/FRA1.chr${chr}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/chr${chr}.sort.BEL1FRA1.badloci
 vcftools --gzvcf ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/chr${chr}.sort.BEL1FRA1.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz

 cat ${outdir}/BEL1.chr${chr}.lmiss ${outdir}/POL3.chr${chr}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/chr${chr}.sort.BEL1POL3.badloci
 vcftools --gzvcf ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/chr${chr}.sort.BEL1POL3.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz

 cat ${outdir}/FRA1.chr${chr}.lmiss ${outdir}/POL3.chr${chr}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/chr${chr}.sort.FRA1POL3.badloci
 vcftools --gzvcf ${outdir}/chr${chr}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/chr${chr}.sort.FRA1POL3.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz

# step 4; print out DP of filtered step 3 and make interval bed file
bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.DP.tsv
Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.female.R  ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.DP.tsv  ${chr} ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed

bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.DP.tsv
Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.female.R  ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.DP.tsv  ${chr} ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed

bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.DP.tsv
Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.female.R  ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.DP.tsv  ${chr} ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed

# step 5, filter AB based on variable sites
#BEL1FRA1
vcftools --gzvcf ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf  | bgzip -c  > ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
tabix ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz 
bcftools concat  --allow-overlaps   ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z -o ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz 
tabix ${outdir}/chr${chr}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

#BEL1POL3
vcftools --gzvcf ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf  | bgzip -c  > ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
tabix ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z -o ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz
tabix ${outdir}/chr${chr}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

#FRA1POL3
vcftools --gzvcf ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf  | bgzip -c  > ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
tabix ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z -o ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz
tabix ${outdir}/chr${chr}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

# # step 6
pair=BEL1FRA1
pixy  --stats pi fst dxy --vcf ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix chr${chr}.${pair}

pair=BEL1POL3
pixy  --stats pi fst dxy --vcf ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix chr${chr}.${pair}

pair=FRA1POL3
pixy  --stats pi fst dxy --vcf ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/chr${chr}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix chr${chr}.${pair}
