#!/bin/bash -l 
#SBATCH -A snic2022-22-105 
#SBATCH -p core  
#SBATCH -n 2
#SBATCH -t 00:10:00  
#SBATCH -J demultiplex_script --mail-type=FAIL,TIME_LIMIT_80

for i in {1..8}
do
    echo $i
    echo Sample_RW${i}
    (echo '#!/bin/bash -l';
	echo "module load bioinfo-tools";
	echo "module load Stacks/2.62";
	echo "mkdir /crex/proj/uppstore2017241/ddRAD_Silene/1-demultiplex/Sample_RW${i}";
	echo "process_radtags -P -p /crex/proj/uppstore2017241/ddRAD_Silene/raw_data_ddRAD_Silene/range_wide_fastq_raw/Sample_RW${i} -b /crex/proj/uppstore2017241/ddRAD_Silene/Range-wide_Xiaodong/range-wide_barcodes/rw_${i}/barcode_rw${i} -o /crex/proj/uppstore2017241/ddRAD_Silene/1-demultiplex/Sample_RW${i}/ -c -q -r -i gzfastq --inline_index --renz_1 ecoRI --renz_2 taqI" ) | sbatch -A snic2022-22-105 -p core -n 8 -t 08:0:0 -J demutiplex --mail-type=FAIL,TIME_LIMIT_80 

done
