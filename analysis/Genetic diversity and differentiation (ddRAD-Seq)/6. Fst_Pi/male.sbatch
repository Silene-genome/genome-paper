#!/bin/bash -l

#SBATCH -p core -n 3                                                                                                                                          
#SBATCH -t 10:00:00                                                                                                                                            
#SBATCH -A naiss2023-5-113

module load bioinfo-tools
module load bcftools
module load vcftools
module load vcflib
module load pixy
module load plink
module load BEDTools

scaff=$1 # eg. chr1, scaffold_1
indir=/crex/proj/uppstore2017241/ddRAD_Silene/4-freebayes/planA_130623/male_haploid_130629
input=${indir}/${scaff}-raw.gvcf
outdir=/crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/result_new/male
output=${outdir}/${scaff}_norm.gvcf


# step 1 norm
vcfallelicprimitives -kg ${input} > ${output} 

# step 2 filter, minDP of 6, remove indel, keep target individuals
vcftools --vcf ${outdir}/${scaff}_norm.gvcf \
  --minDP 6 \
  --remove-indels \
  --recode-INFO-all \
  --keep ${outdir}/sample_goodDP.male.tsv  \
  --recode --stdout | bgzip -c > ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz

# # step 3 filter data missingness per pop
  vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --keep ${outdir}/BEL1_goodDP.list --missing-site --out ${outdir}/BEL1.${scaff}
  vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --keep ${outdir}/FRA1_goodDP.list --missing-site --out ${outdir}/FRA1.${scaff}
  vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --keep ${outdir}/POL3_goodDP.list --missing-site --out ${outdir}/POL3.${scaff}

 thres=0.5
 cat ${outdir}/BEL1.${scaff}.lmiss ${outdir}/FRA1.${scaff}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/${scaff}.sort.BEL1FRA1.badloci
 vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/${scaff}.sort.BEL1FRA1.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz

 cat ${outdir}/BEL1.${scaff}.lmiss ${outdir}/POL3.${scaff}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/${scaff}.sort.BEL1POL3.badloci
 vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/${scaff}.sort.BEL1POL3.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz

 cat ${outdir}/FRA1.${scaff}.lmiss ${outdir}/POL3.${scaff}.lmiss | awk '!/CHR/' | awk '$6 > 0.5' | sort -k2 -n | uniq > ${outdir}/${scaff}.sort.FRA1POL3.badloci
 vcftools --gzvcf ${outdir}/${scaff}_rmlowDPsample_noINDEL_minDP6.gvcf.gz  --exclude-positions ${outdir}/${scaff}.sort.FRA1POL3.badloci  --recode --recode-INFO-all --stdout | bgzip -c > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz


# step 4; print out DP of filtered step 3 and make interval bed file
  bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.DP.tsv
  Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.male.R  ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.DP.tsv  ${scaff} ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed

  bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.DP.tsv
  Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.male.R  ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.DP.tsv  ${scaff} ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed

  bcftools query -f '%CHROM\t%POS[\t%DP]\n' ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz   > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.DP.tsv
  Rscript /crex/proj/uppstore2017241/ddRAD_Silene/5-vcffilter/script_new/makeInterval.male.R  ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.DP.tsv  ${scaff} ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed


# step 5, filter Allelic Balance based on variable sites
# #BEL1FRA1          
vcftools --gzvcf ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf
tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z -o ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz
tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

# #BEL1POL3
vcftools --gzvcf ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf  | bgzip -c  > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z -o ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz
tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

# #FRA1POL3 
vcftools --gzvcf ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --max-maf 0  --recode-INFO-all  --recode --stdout | bgzip -c > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz
vcftools --gzvcf ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.gvcf.gz --mac 1  --recode-INFO-all  --recode --stdout  > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf  | bgzip -c  > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  && tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.variant.vcf.gz  -O z  -o ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz
tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz

# # step 6


 pair=BEL1FRA1
 pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}

 pair=BEL1POL3
 pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}

 pair=FRA1POL3
 pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_AB25.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}


#step 5 and 6 for Y scaff

#BEL1FRA1
bgzip -c ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf > ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf.gz
tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf.gz
plink  --vcf ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf --allow-extra-chr --make-bed --out ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant
plink --bfile ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant  --allow-extra-chr   --set-missing-var-ids "@:#:\$1:\$2"  --het --hardy --ibc --missing --out ${outdir}/${scaff}_BEL1FRA1
awk '$7<=0.1 {print $2}' ${outdir}/${scaff}_BEL1FRA1.hwe | sed '1d' > ${outdir}/${scaff}_BEL1FRA1.OHE01.keep.sites
cut ${outdir}/${scaff}_BEL1FRA1.OHE01.keep.sites -d":" -f1,2 | sed 's/:/\t/g' > ${outdir}/${scaff}_BEL1FRA1.OHE01.keep.sites.regionfile
bcftools view --regions-file ${scaff}_BEL1FRA1.OHE01.keep.sites.regionfile ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.variant.vcf.gz  -Oz -o ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz  -O z -o ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz 
tabix ${outdir}/${scaff}_BEL1FRA1_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz

pair=BEL1FRA1  
pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}


#BEL1POL3
bgzip -c ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf > ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz
tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz
plink  --vcf ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf --allow-extra-chr --make-bed --out ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant
plink --bfile ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant  --allow-extra-chr   --set-missing-var-ids "@:#:\$1:\$2"  --het --hardy --ibc --missing --out ${outdir}/${scaff}_BEL1POL3
awk '$7<=0.1 {print $2}' ${outdir}/${scaff}_BEL1POL3.hwe | sed '1d' > ${outdir}/${scaff}_BEL1POL3.OHE01.keep.sites
cut ${outdir}/${scaff}_BEL1POL3.OHE01.keep.sites -d":" -f1,2 | sed 's/:/\t/g' > ${outdir}/${scaff}_BEL1POL3.OHE01.keep.sites.regionfile
bcftools view --regions-file ${scaff}_BEL1POL3.OHE01.keep.sites.regionfile ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz  -Oz -o ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz  -O z -o ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz 
tabix ${outdir}/${scaff}_BEL1POL3_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz

pair=BEL1POL3  
pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}

#FRA1POL3
bgzip -c ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf > ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz
tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz
plink  --vcf ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf --allow-extra-chr --make-bed --out ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant
plink --bfile ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant  --allow-extra-chr   --set-missing-var-ids "@:#:\$1:\$2"  --het --hardy --ibc --missing --out ${outdir}/${scaff}_FRA1POL3
awk '$7<=0.1 {print $2}' ${outdir}/${scaff}_FRA1POL3.hwe | sed '1d' > ${outdir}/${scaff}_FRA1POL3.OHE01.keep.sites
cut ${outdir}/${scaff}_FRA1POL3.OHE01.keep.sites -d":" -f1,2 | sed 's/:/\t/g' > ${outdir}/${scaff}_FRA1POL3.OHE01.keep.sites.regionfile
bcftools view --regions-file ${scaff}_FRA1POL3.OHE01.keep.sites.regionfile ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.variant.vcf.gz  -Oz -o ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz
bcftools concat  --allow-overlaps   ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop.invariant.vcf.gz  ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_OHE01.variant.vcf.gz  -O z -o ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz 
tabix ${outdir}/${scaff}_FRA1POL3_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz

pair=FRA1POL3  
pixy  --stats pi fst dxy --vcf ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop_OHE01.allsites.vcf.gz  --populations ${outdir}/sample_goodDP.${pair}.popfile --bed_file  ${outdir}/${scaff}_${pair}_noINDEL_minDP6_geno05perPop.minInteral100.pixy.bed --fst_type hudson  --output_folder ${outdir} --output_prefix ${scaff}.${pair}


