#!/bin/bash


# ./downsample_data.sh DS_RATIO INPUT_DIR PREFIX SUFFIX OUTPUT_DIR

#./downsample_data.sh 16 /mnt/pd2/shared/EZH2_Interactome/KARPAS/DMSO/Separated_filtered KARPAS_DMSO .rmdup.threshold0_mapq1.txt OUTPUT_DOWNSAMPLE

# frag_file <- file.path(setDir, "/mnt/pd2/shared/EZH2_Interactome/KARPAS/DMSO/Separated_filtered/KARPAS_DMSO_chr1_library1_Rep1AB.rmdup.threshold0_mapq1.txt")


start_time=$(date -R)    
set -e

if [ "$#" -ne 5 ]; then
    echo "Illegal number of parameters"
    exit
fi



ds_ratio="$1"
in_dir="$2"
file_prefix="$3"
file_suffix="$4"
out_dir="$5"

all_chrs=( {1..22} )

echo "********************** START"
echo "in_dir = $in_dir"
echo "out_dir = $out_dir/$file_prefix"

output_dir="$out_dir/$file_prefix"
mkdir -p $output_dir


for chrom in ${all_chrs[@]}; do

    chr="chr${chrom}"

    logFile="$output_dir/${file_prefix}_${chr}_downsample_logFile.txt"
    rm -f logFile

    echo "> START $chr"
    echo "... logFile = $logFile"

    #all_files=$(( ls $in_dir/${file_prefix}_${chr}_*${file_suffix} ))
    all_files=`ls $in_dir/${file_prefix}_${chr}_*${file_suffix}`

    echo "Found fragment files:" >> $logFile
    wc -l $in_dir/${file_prefix}_${chr}_*${file_suffix} >> $logFile

    # DOWNSAMPLE
    for frag_file in ${all_files[@]}; do
        echo "> start $frag_file"
        num=$(cat $frag_file | wc -l)
        num_downsample=`expr $(($num/$ds_ratio))`
        echo "... downsample from $num -> $num_downsample"
        out_name=`basename $frag_file $file_suffix`
        out_file="$output_dir/${out_name}_downsample${ds_ratio}.txt"
        shuf -n $num_downsample $frag_file > $out_file
        echo "... written: $out_file"
    done

    echo "Written downsampled files:" >> $logFile
    wc -l $output_dir/*_downsample${ds_ratio}.txt >> $logFile

    # MERGE ALL DOWNSAMPLED FILES (AND SORT)
    out_merged_file="$out_dir/${file_prefix}_${chr}_downsample${ds_ratio}_merged.txt"
    cat $output_dir/${file_prefix}_${chr}*downsample${ds_ratio}.txt | sort -k3,3n -k6,6n > $out_merged_file
    echo "... written: $out_merged_file"
    echo "Written merged file with DS:" >> $logFile
    wc -l $out_merged_file >> $logFile

    # MERGE THE NOT DOWNSAMPLED FILES (AND SORT)
    out_merged_file_noDS="$out_dir/${file_prefix}_${chr}_noDS_merged.txt"
    cat  $in_dir/${file_prefix}_${chr}_*${file_suffix}| sort -k3,3n -k6,6n > $out_merged_file_noDS
    echo "... written: $out_merged_file_noDS"
    echo "Written merged file no DS:" >> $logFile
    wc -l $out_merged_file_noDS >> $logFile


    echo "written: $logFile"

done




########## END ####################################################################################################################################


echo "*** DONE"
end_time=$(date -R)    
echo $start_time
echo $end_time

exit 0

# PPATH=$(dirname $(readlink -f "$0")) # to retrieve directory from where the script is launched


#resolution=$2
#juicer_path=$3
#mkdir -p "$DPATH/$CELL"
##merge raw data 
#find "$DPATH/$CELL" -name "*_merged_nodups.txt.gz"|xargs zcat | sort -k3,3d -k7,7d > "$DPATH/$CELL/total_merged_nodups.txt"
##downsample
#ratio=16
#num=$(cat $DPATH/$CELL/total_merged_nodups.txt |wc -l)
#num_downsample=`expr $(($num/$ratio))`
#shuf -n $num_downsample $DPATH/$CELL/total_merged_nodups.txt | sort -k3,3d -k7,7d  > $DPATH/$CELL/total_merged_nodups_downsample_ratio_$ratio.txt
#echo "merge data done!"

#mkdir -p "$DPATH/$CELL/intra_NONE"
#mkdir -p "$DPATH/$CELL/intra_VC"
#mkdir -p "$DPATH/$CELL/intra_KR"

##write your own path to juicer tools here
#juicer_tool="/home/liuqiao/software/juicer/scripts/juicer_tools.jar"
##generate .HIC file using juicer tool, -xmx50g indicates 50g for memory which can be replaced with an appropriate value
#java -Xmx50g  -jar $juicer_tool pre $DPATH/$CELL/total_merged_nodups.txt $DPATH/$CELL/total_merged.hic hg19
#java -Xmx50g  -jar $juicer_tool pre $DPATH/$CELL/total_merged_nodups_downsample_ratio_$ratio.txt $DPATH/$CELL/total_merged_downsample_ratio_$ratio.hic hg19

#chromes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 "X" "Y")

##generate Hi-C raw contacts using iuicer tool
#for chrom in ${chromes[@]}
#do
#	bash $PPATH/sbatch_juicer_script.sh $chrom $DPATH $CELL $resolution $juicer_path
#	#replace with "sbatch sbatch_juicer_script.sh $chrom $DPATH $CELL $resolution $juicer_path" if slurm is installed, this will save a lot time.
#done
