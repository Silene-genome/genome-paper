#!/bin/bash 

source ../config/config
source ../config/colors

echo rnaseqlist is $RNAseqlist 
#microscript to run all rnaseq steps from read trimming to read mapping and count
#RNAseq=YES

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="log_files"

#create folder if not existent:
mkdir $LOG_FOLDER 2>/dev/null

haplotype=$1

genome=03_genome/"$haplotype".fa

#check that the file in 03_genome exist:

#check that the raw_data for RNAseq are present:

if [ -f Trimmomatic-0.39/trimmomatic-0.39.jar ] ; then 
    echo "found trimmomatic jar " ; 
else 
    echo -e "trimmomatic jar not found\nwill attempt to download" 
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip 
    if ! unzip Trimmomatic-0.39.zip
    then
        echo -e "${RED} ERROR! failed downloading trimmomatic - check your internet connexion  \n${NC}"
        exit 1
    else
        echo -e "trimmomatic download successfull\n"
fi

fi

#launch gmap :
../00_scripts/02_gmap.sh "$genome" 

echo -e "\n\ntrimming read for RNAseq\n" 
#paste <(ls 01_raw_data/*1.f**q.gz ) <(ls 01_raw_data/*2.f**q.gz ) > file1file2.tmp

ncol=$(awk '{print NF}'  "$RNAseqlist" |uniq)

if [[ $ncol = 2 ]] ; then
    #assuming PE:
    echo -e "\nrunning trimmomatic assuming reads are Paired-End\n" 
    while IFS=$'\t' read -r -a read ; 
    do 
        ../00_scripts/01_trimmomatic_PE.sh "${read[0]}" "${read[1]}"  
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR TRIMMOMATIC FAILED"
            exit 1
        fi

    done < "$RNAseqlist" #  file1file2.tmp 

    if [ $? -eq 0 ]; then
        echo "trimmomatic complete "
        echo -e "\ncounting the number of retained reads\n"        #rm file1file2.tmp 
        ../00_scripts/utility_scripts/count_read_fastq.sh 02_trimmed/*gz > read_count.txt

    else
        echo -e "\n${RED}#ERROR : Runnning trimmomatic failed. please check your input files${NC}"
        exit 1
    fi

    #launch gsnap - samtools and read count:
    for read1 in 02_trimmed/*R1.paired.fastq.gz  ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_PE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi

    done 

else
    #assuming SE:
    echo "running trimmomatic" 
    while IFS=$'\t' read -r -a read ; 
    do 
        ../00_scripts/01_trimmomatic_SE.sh "${read[0]}" 2>&1 |tee "$LOG_FOLDER"/trimmo_"${read[0]}"_log  
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR TRIMMOMATIC FAILED"
            exit 1
        fi

    done < "$RNAseqlist" #  file1file2.tmp 
    
    if [ $? -eq 0 ]; then
        echo trimmomatic complete
        #rm file1file2.tmp 
    else
        echo -e "\n${RED} ERROR : Runnning trimmomatic failed. please check your input files ${NC}"
        exit 1
    fi

    #launch gsnap - samtools and read count:
    for read1 in 02_trimmed/*R1.paired.fastq.gz ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_SE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi

    done 

fi

rm -rf Trimmomatic-0.39* 
