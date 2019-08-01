#!/bin/bash


# ./all_prep_hicgan.sh  # no args, everything hard-coded

start_time=$(date -R)    
set -e

step1=1
step2=1
step3=1

bin_size=10000
rexec="Rscript"
pyexec="python"

# STEP1: downsample and merge count files -> purely bash
ds_ratio=16
in_dir_step1="/mnt/pd2/shared/EZH2_Interactome/KARPAS/DMSO/Separated_filtered"
file_prefix_step1="KARPAS_DMSO"
out_dir_step1="INPUT_FRAGS"
frag_suffix=".rmdup.threshold0_mapq1.txt"        
output_dir_step1="$out_dir_step1/$file_prefix_step1"

# STEP2: aggregate the counts
step2_script="2_agg_frags.R"
out_dir_step2="INPUT_AGG"
output_dir_step2="$out_dir_step2/$file_prefix_step1"

# STEP3: split and prepare data for hicGAN
# python data_split.py <INPUT_DIR> <FILE_PREFIX> <FILE_SUFFIX HIGH-RESOL> <FILE_SUFFIX LOW-RESOL> <OUT_DIR>
# python data_split.py INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt downsample16_merged_agg.txt INPUT_MATS
step3_script="data_split.py"
out_dir_step3="INPUT_HKL"
hrSuffix_step3="noDS_merged_agg.txt"
lrSuffix_step3="downsample16_merged_agg.txt"
output_dir_step3="$out_dir_step3/$file_prefix_step1"



echo "!!! HARD-CODED:"
echo "... bin_size = $bin_size"
echo "... ds_ratio = $ds_ratio"


all_chrs=( {1..22} )


for chrom in ${all_chrs[@]}; do

    chr="chr${chrom}"

    #*************************************************
    #**** STEP1: downsample and merge count files
    #*************************************************
    if [[ $step1 -eq 1 ]] ; then



        mkdir -p $output_dir_step1


        logFile_step1="$output_dir_step1/${file_prefix_step1}_${chr}_downsample_logFile.txt"
        rm -f $logFile_step1

        echo "> START $chr - STEP1"
        echo "... in_dir_step1_step1 = $in_dir_step1_step1"
        echo "... output_dir_step1 = $output_dir_step1"
        echo "... logFile_step1 = $logFile_step1"


        #all_files=$(( ls $in_dir_step1/${file_prefix_step1}_${chr}_*${frag_suffix} ))
        all_files=`ls $in_dir_step1/${file_prefix_step1}_${chr}_*${frag_suffix}`

        echo "Found fragment files:" >> $logFile_step1
        wc -l $in_dir_step1/${file_prefix_step1}_${chr}_*${frag_suffix} >> $logFile_step1

        # DOWNSAMPLE
        for frag_file in ${all_files[@]}; do
            echo "> start $frag_file"
            num=$(cat $frag_file | wc -l)
            num_downsample=`expr $(($num/$ds_ratio))`
            echo "... downsample from $num -> $num_downsample"
            out_name=`basename $frag_file $frag_suffix`
            out_file="$output_dir_step1/${out_name}_downsample${ds_ratio}.txt"
            shuf -n $num_downsample $frag_file > $out_file
            echo "... written: $out_file"
        done

        echo "Written downsampled files:" >> $logFile_step1
        wc -l $output_dir_step1/*_downsample${ds_ratio}.txt >> $logFile_step1

        # MERGE ALL DOWNSAMPLED FILES (AND SORT)
        cat "... merge downsampled data"
        out_merged_file="$output_dir_step1/${file_prefix_step1}_${chr}_downsample${ds_ratio}_merged.txt"
        cat $output_dir_step1/${file_prefix_step1}_${chr}*downsample${ds_ratio}.txt | sort -k3,3n -k6,6n > $out_merged_file
        echo "... written: $out_merged_file"
        echo "Written merged file with DS:" >> $logFile_step1
        wc -l $out_merged_file >> $logFile_step1

        # MERGE THE NOT DOWNSAMPLED FILES (AND SORT)
        cat "... merge init data"
        out_merged_file_noDS="$output_dir_step1/${file_prefix_step1}_${chr}_noDS_merged.txt"
        cat  $in_dir_step1/${file_prefix_step1}_${chr}_*${frag_suffix}| sort -k3,3n -k6,6n > $out_merged_file_noDS
        echo "... written: $out_merged_file_noDS"
        echo "Written merged file no DS:" >> $logFile_step1
        wc -l $out_merged_file_noDS >> $logFile_step1


        echo "written: $logFile_step1"

        # OUTPUT FILES LOOK LIKE:
        # INPUT_FRAGS/KARPAS_DMSO/KARPAS_DMSO_chr1_downsample16_merged.txt
        # INPUT_FRAGS/KARPAS_DMSO/KARPAS_DMSO_chr1_noDS_merged.txt

    fi # end-if STEP1


    #*************************************************
    #**** STEP2: aggregate the counts 
    #*************************************************
    if [[ $step2 -eq 1 ]] ; then

        echo "> START $chr - STEP2"



        mkdir -p $output_dir_step2


        logFile_step2="$output_dir_step2/${file_prefix_step1}_${chr}_aggreg_logFile.txt"
        echo "... logFile_step2 = $logFile_step2"
        rm -f $logFile_step2


        step1_dsFile="$output_dir_step1/${file_prefix_step1}_${chr}_downsample${ds_ratio}_merged.txt"
        ds_cmd="$rexec $step2_script $step1_dsFile $bin_size $output_dir_step2"
        
        echo "> $ds_cmd" >> $logFile_step2
        $ds_cmd >> $logFile_step2

        step1_nodsFile="$output_dir_step1/${file_prefix_step1}_${chr}_noDS_merged.txt"
        nods_cmd="$rexec $step2_script $step1_nodsFile $bin_size $output_dir_step2"
        echo "> $nods_cmd" >> $logFile_step2
        $nods_cmd >> $logFile_step2

        # output of step2 is then (hard-coded stuff in Rscript): 
        # "$output_dir_step2/${file_prefix_step1}_${chr}_noDS_merged_agg.txt"
        # "$output_dir_step2/${file_prefix_step1}_${chr}_downsample${ds_ratio}_merged_agg.txt"

        echo "written: $logFile_step2"


        # OUTPUT FILES LOOK LIKE:
        # INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_downsample16_merged_agg.txt
        # INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_noDS_merged_agg.txt


    fi # end-if STEP2




done # end-iterating over chromo



#*************************************************
#**** STEP3: prepare matrix for hicGAN
#*************************************************
if [[ $step3 -eq 1 ]] ; then

    echo "> START - STEP3"

    mkdir -p $output_dir_step3

    logFile_step3="$output_dir_step3/${file_prefix_step1}_prepData_logFile.txt"
    echo "... logFile_step3 = $logFile_step3"
    rm -f $logFile_step3


    # python data_split.py INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt downsample16_merged_agg.txt INPUT_MATS
    $pyexec $step3_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 $lrSuffix_step3 $output_dir_step3 >> $logFile_step3

    echo "written: $logFile_step3"

fi # end-if STEP3



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
