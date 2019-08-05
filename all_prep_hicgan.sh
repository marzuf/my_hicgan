#!/bin/bash


# ./all_prep_hicgan.sh  # no args, everything hard-coded

start_time=$(date -R)    
#set -e

step1=0 # downsample and merge count files (bash)
step2=0 # aggregate counts (R script)
step3=0 # prepare data for hicGAN (python script)
step4=0 # run hicGAN (python script)
step5a=0 # prepare data for FitHiC, low-resol and high-resol data (R script)
step5b=0 # prepare data for FitHiC, super-resol data (R script)
step6a=1 # run FitHiC for low-resol and high-resol data (R script)
step6b=0 # run FitHiC for super-resol data (R script)

bin_size=10000
bin_sizeKb=$((bin_size / 1000))

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
step3_script="3_data_split.py"
out_dir_step3="INPUT_HKL"
hrSuffix_step3="${bin_sizeKb}kb_noDS_merged_agg.txt"
lrSuffix_step3="${bin_sizeKb}kb_downsample${ds_ratio}_merged_agg.txt"
output_dir_step3="$out_dir_step3/$file_prefix_step1"

# STEP5: prepare files for FitHiC
step5_script="5_prep_FitHiC.R"
out_dir_step5="PREP_FITHIC"
output_dir_step5="$out_dir_step5/$file_prefix_step1"

# STEP4a: train hicGAN
step4a_script="4_train_hicGAN"
# STEP4b: evaluate hicGAN
# STEP4c: build super-resol data

# STEP6: run FitHiC
step6_script="6_run_FitHiC.R"
out_dir_step6="RUN_FITHIC"
output_dir_step6="$out_dir_step6/$file_prefix_step1"
hrSuffix_step5="${bin_sizeKb}kb_noDS_merged_agg"
lrSuffix_step5="${bin_sizeKb}kb_downsample${ds_ratio}_merged_agg"


echo "!!! HARD-CODED:"
echo "... bin_size = $bin_size"
echo "... ds_ratio = $ds_ratio"


all_chrs=( {1..22} )
#all_chrs=( 1 )


for chrom in ${all_chrs[@]}; do

    chr="chr${chrom}"

    #*************************************************
    #**** STEP1: downsample and merge count files
    #*************************************************
    if [[ $step1 -eq 1 ]] ; then



        mkdir -p $output_dir_step1


        logFile_step1="$output_dir_step1/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_downsample_logFile.txt"
        rm -f $logFile_step1

        echo "> START $chr - STEP1"
        echo "... in_dir_step1 = $in_dir_step1"
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
            out_file="$output_dir_step1/${out_name}_${bin_sizeKb}kb_downsample${ds_ratio}.txt"
            shuf -n $num_downsample $frag_file > $out_file
            echo "... written: $out_file"
        done

        echo "Written downsampled files:" >> $logFile_step1
        wc -l $output_dir_step1/*_downsample${ds_ratio}.txt >> $logFile_step1

        # MERGE ALL DOWNSAMPLED FILES (AND SORT)
        echo "... merge downsampled data"
        out_merged_file="$output_dir_step1/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_downsample${ds_ratio}_merged.txt"
        cat $output_dir_step1/${file_prefix_step1}_${chr}*downsample${ds_ratio}.txt | sort -k3,3n -k6,6n > $out_merged_file
        echo "... written: $out_merged_file"
        echo "Written merged file with DS:" >> $logFile_step1
        wc -l $out_merged_file >> $logFile_step1

        # MERGE THE NOT DOWNSAMPLED FILES (AND SORT)
        echo "... merge init data"
        out_merged_file_noDS="$output_dir_step1/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_noDS_merged.txt"
        cat  $in_dir_step1/${file_prefix_step1}_${chr}_*${frag_suffix}| sort -k3,3n -k6,6n > $out_merged_file_noDS
        echo "... written: $out_merged_file_noDS"
        echo "Written merged file no DS:" >> $logFile_step1
        wc -l $out_merged_file_noDS >> $logFile_step1


        echo "written: $logFile_step1"

        # OUTPUT FILES LOOK LIKE:
        # INPUT_FRAGS/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged.txt
        # INPUT_FRAGS/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged.txt

    fi # end-if STEP1


    #*************************************************
    #**** STEP2: aggregate the counts 
    #*************************************************
    if [[ $step2 -eq 1 ]] ; then

        echo "> START $chr - STEP2"



        mkdir -p $output_dir_step2


        logFile_step2="$output_dir_step2/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_aggreg_logFile.txt"
        echo "... logFile_step2 = $logFile_step2"
        rm -f $logFile_step2


        step1_dsFile="$output_dir_step1/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_downsample${ds_ratio}_merged.txt"
        ds_cmd="$rexec $step2_script $step1_dsFile $bin_size $output_dir_step2"
        
        echo "> $ds_cmd" >> $logFile_step2
        $ds_cmd >> $logFile_step2

        step1_nodsFile="$output_dir_step1/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_noDS_merged.txt"
        nods_cmd="$rexec $step2_script $step1_nodsFile $bin_size $output_dir_step2"
        echo "> $nods_cmd" >> $logFile_step2
        $nods_cmd >> $logFile_step2

        # output of step2 is then (hard-coded stuff in Rscript): 
        # "$output_dir_step2/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_noDS_merged_agg.txt"
        # "$output_dir_step2/${file_prefix_step1}_${chr}_${bin_sizeKb}kb_downsample${ds_ratio}_merged_agg.txt"

        echo "written: $logFile_step2"


        # OUTPUT FILES LOOK LIKE:
        # INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged_agg.txt
        # INPUT_AGG/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged_agg.txt


    fi # end-if STEP2




done # end-iterating over chromo



#*************************************************
#**** STEP3: prepare matrix for hicGAN
#*************************************************
if [[ $step3 -eq 1 ]] ; then

    echo "> START - STEP3"

    mkdir -p $output_dir_step3

    logFile_step3="$output_dir_step3/${file_prefix_step1}_prepData_hicGAN_logFile.txt"
    echo "... logFile_step3 = $logFile_step3"
    rm -f $logFile_step3


    # python data_split.py INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt downsample16_merged_agg.txt INPUT_HKL
    echo $pyexec $step3_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 $lrSuffix_step3 $output_dir_step3 
    $pyexec $step3_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 $lrSuffix_step3 $output_dir_step3 >> $logFile_step3

    echo "written: $logFile_step3"

fi # end-if STEP3



#*************************************************
#**** STEP4a: train hicGAN (train hicGAN)
#*************************************************

if [[ $step4a -eq 1 ]] ; then

    echo ""


fi # end-if STEP4a



#*************************************************
#**** STEP4b: evaluate hicGAN (evaluate hicGAN)
#*************************************************

if [[ $step4b -eq 1 ]] ; then

    echo ""
fi # end-if STEP4b


#*************************************************
#**** STEP4c: run hicGAN (build super-resol matrix)
#*************************************************

if [[ $step4c -eq 1 ]] ; then

    echo ""


fi # end-if STEP4c


#*************************************************
#**** STEP5a: prepare data for FitHiC (low-resol and high-resol)
#*************************************************

# Rscript 5_prep_FitHiC.R INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt chr1 10000 PREP_FITHIC

if [[ $step5a -eq 1 ]] ; then


    for chrom in ${all_chrs[@]}; do

        chr="chr${chrom}"


        echo "> START - STEP5a"

        mkdir -p $output_dir_step5

        logFile_step5a="$output_dir_step5/${file_prefix_step1}_prepData_FitHiC_logFile.txt"
        echo "... logFile_step5a = $logFile_step5a"
        rm -f $logFile_step5a

   
        # Rscript 5_prep_FitHiC.R INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt chr1 10000 PREP_FITHIC/KARPAS_DMSO
        echo "$rexec $step5_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 $chr $bin_size $output_dir_step5"
        $rexec $step5_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 $chr $bin_size $output_dir_step5 >> $logFile_step5a

        echo "$rexec $step5_script $output_dir_step2 $file_prefix_step1 $lrSuffix_step3 $chr $bin_size $output_dir_step5"
        $rexec $step5_script $output_dir_step2 $file_prefix_step1 $lrSuffix_step3 $chr $bin_size $output_dir_step5 >> $logFile_step5a





        echo "written: $logFile_step5a"

        # OUTPUT FILES LOOK LIKE:
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged_agg_FitHiC_fragsfile.txt
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged_agg_FitHiC_intersfile.txt
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged_agg_FitHiC_fragsfile.txt
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_downsample16_merged_agg_FitHiC_intersfile.txt


    done


fi # end-if STEP5a


#*************************************************
#**** STEP5b: prepare data for FitHiC (super-resol)
#*************************************************
if [[ $step5b -eq 1 ]] ; then


    for chrom in ${all_chrs[@]}; do

        chr="chr${chrom}"


        echo "> START - STEP5b"

        mkdir -p $output_dir_step5

        logFile_step5b="$output_dir_step5/${file_prefix_step1}_prepData_superResol_FitHiC_logFile.txt"
        echo "... logFile_step5b = $logFile_step5b"
        rm -f $logFile_step5b

   
        # Rscript 5_prep_FitHiC.R INPUT_AGG/KARPAS_DMSO KARPAS_DMSO noDS_merged_agg.txt chr1 10000 PREP_FITHIC/KARPAS_DMSO
        #echo "$rexec $step5_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 chr$chrom $bin_size $output_dir_step5"
        #$rexec $step5_script $output_dir_step2 $file_prefix_step1 $hrSuffix_step3 chr$chrom $bin_size $output_dir_step5 >> $logFile_step5a




        echo "written: $logFile_step5b"

        # OUTPUT FILES LOOK LIKE:
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_XXX_FitHiC_fragsfile.txt
        # PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_XXX_FitHiC_intersfile.txt


    done


fi # end-if STEP5b


#*************************************************
#**** STEP6a: run FitHiC (low-resol and high-resol)
#*************************************************
if [[ $step6a -eq 1 ]] ; then



    for chrom in ${all_chrs[@]}; do

        chr="chr${chrom}"


        echo "> START - STEP6a"

        mkdir -p $output_dir_step6

        logFile_step6a="$output_dir_step6/${file_prefix_step1}_prepData_runFitHiC_logFile.txt"
        echo "... logFile_step6a = $logFile_step6a"
        rm -f $logFile_step6a


        # Rscript 6_run_FitHiC.R PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged_agg_FitHiC_fragsfile.txt PREP_FITHIC/KARPAS_DMSO/KARPAS_DMSO_chr1_10kb_noDS_merged_agg_FitHiC_intersfile.txt RUN_FITHIC/KARPAS_DMSO/low_resol

        # run for high resol data:
        hr_fragfile="${output_dir_step5}/${file_prefix_step1}_${chr}_${hrSuffix_step5}_FitHiC_fragsfile.txt"
        hr_interfile="${output_dir_step5}/${file_prefix_step1}_${chr}_${hrSuffix_step5}_FitHiC_intersfile.txt"
        echo $rexec $step6_script $hr_fragfile $hr_interfile $output_dir_step6/high_resol/$chr
        $rexec $step6_script $hr_fragfile $hr_interfile $output_dir_step6/high_resol/$chr >> $logFile_step6a


        # run for low resol data:
        lr_fragfile="${output_dir_step5}/${file_prefix_step1}_${chr}_${lrSuffix_step5}_FitHiC_fragsfile.txt"
        lr_interfile="${output_dir_step5}/${file_prefix_step1}_${chr}_${lrSuffix_step5}_FitHiC_intersfile.txt"
        echo $rexec $step6_script $lr_fragfile $lr_interfile $output_dir_step6/low_resol/$chr
        $rexec $step6_script $lr_fragfile $lr_interfile $output_dir_step6/low_resol/$chr >> $logFile_step6a

        echo "written: $logFile_step6a"

        # OUTPUT FILES LOOK LIKE:

    done





fi # end-if STEP6a


#*************************************************
#**** STEP6b: run FitHiC (super-resol)
#*************************************************
if [[ $step6b -eq 1 ]] ; then



    for chrom in ${all_chrs[@]}; do

        chr="chr${chrom}"


        echo "> START - STEP6b"

        mkdir -p $output_dir_step6

        logFile_step6b="$output_dir_step6/${file_prefix_step1}_prepData_superResol_runFitHiC_logFile.txt"
        echo "... logFile_step6b = $logFile_step6b"
        rm -f $logFile_step6b


        echo "written: $logFile_step6b"


    done


fi # end-if STEP6b






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
