#!/bin/bash

# Purpose: Collect stats from FastQC, Trimmomatic, STAR into one table
# For RNAseq data
# Currently for PE reads
# Has relative paths and assumes directory structure:
#   -ECP_Project_dir
#     - scripts
#     - data
#       - FastQC
#           - raw
#       - star
#     - slurm_debug
# execute script in scripts directory and final output is in ECP_Project_dir    

#-----------------------------------------------------------------------

# Current directory is scripts directory
dir=$(pwd)
# change to project home dir
cd ..

#-----------#
# Variables #
#-----------#
# path to raw fastq files (which was used for trimmomatic)
# goal is to remove this string from the trimmomatic log, leaving only fastq file name
raw_file_path=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP59/data
# suffix of fastq files to remove
suffix=_R1.fastq.gz
# suffix of unzipped sample fastqc folder to remove, leaving fastq prefix for merging with metadata
#   for now, manually edit line 69 so sed statement is correct. Usually "\_R1\_001\_fastqc\" or "\_R1\_fastqc\" etc
#   includes the fastq file suffix (excluding the fastq.gz) and the name of the fastqc_data.txt file and "Total Sequences"
#fqc_folder_suffix=
# path to dir containing FastQC raw results
fqc_path=./data/FastQC/raw
# path to dir containing trimomatic logs (typically in SLURM .err files)
trim_path=./slurm_debug
# path to dir containing STAR logs
star_path=./data/star
# path to metadata xlsx file
metadata=./metadata/ECP59.metdata.v2.xlsx

#------------------------------------------------------------------------

#-----------------------#
# Get Stats from FastQC #
#-----------------------#
cd $fqc_path
#for f in `ls --color=none *.zip`; do unzip $f; done
# Grab Metrics
grep "Total Sequences" */fastqc_data.txt > tmp_reads.txt
grep "Sequence length" */fastqc_data.txt > tmp_readLength.txt
grep "Total Deduplicated" */fastqc_data.txt > tmp_duplicated.txt
grep "^%GC" */fastqc_data.txt > tmp_GC.txt
# Clean up Metrics
### gather the R2 entries into a separate file
awk 'NR%2==1' tmp_reads.txt > tmp_readsR1.txt
awk 'NR%2==0' tmp_reads.txt > tmp_readsR2.txt
awk 'NR%2==1' tmp_readLength.txt > tmp_readLengthR1.txt
awk 'NR%2==0' tmp_readLength.txt > tmp_readLengthR2.txt
awk 'NR%2==1' tmp_GC.txt > tmp_GCR1.txt
awk 'NR%2==0' tmp_GC.txt > tmp_GCR2.txt
### transform the deduplicated percentage to duplicated percentage (by subtracting from 100)
awk -F'\t' '{print $1,(100-$2)}' tmp_duplicated.txt > tmp_duplicated2.txt
### gather the R2 entries into a separate file
awk 'NR%2==1' tmp_duplicated2.txt > tmp_duplicatedR1.txt
awk 'NR%2==0' tmp_duplicated2.txt > tmp_duplicatedR2.txt
### Keep only the fastq prefix name in the first field of the tmp_readsR1.txt file
sed -i 's/\_R1\_fastqc\/fastqc_data\.txt\:Total Sequences/''/g' tmp_readsR1.txt
### Keep only the last field for remaining tmp R1,R2 files
cut -f2 tmp_readsR2.txt > tmp2_readsR2.txt
cut -f2 tmp_readLengthR1.txt > tmp2_readLengthR1.txt
cut -f2 tmp_readLengthR2.txt > tmp2_readLengthR2.txt
cut -f2 tmp_GCR1.txt > tmp2_GCR1.txt
cut -f2 tmp_GCR2.txt > tmp2_GCR2.txt
cut -d ' ' -f4 tmp_duplicatedR1.txt > tmp2_duplicatedR1.txt
cut -d ' ' -f4 tmp_duplicatedR2.txt > tmp2_duplicatedR2.txt
# Combine FQC stats into one table
paste -d "\t" tmp_readsR1.txt tmp2_readsR2.txt tmp2_readLengthR1.txt tmp2_readLengthR2.txt tmp2_GCR1.txt tmp2_GCR2.txt tmp2_duplicatedR1.txt tmp2_duplicatedR2.txt > tmp_fqc_stats.tmp.txt
sort tmp_fqc_stats.tmp.txt > fqc_stats.table.txt

# add a header
sed -i '1s/^/Fastq_File\tR1_raw_reads\tR2_raw_reads\tR1_read_length\tR2_read_length\tR1_GC%\tR2_GC%\tR1_dup%\tR2_dup%\n/' fqc_stats.table.txt
rm tmp*.txt
# Move output file and change directories
cd $dir
cd ..
mv $fqc_path/fqc_stats.table.txt .

#----------------------------------------------------------------------------------#

#----------------------------#
# Get Stats From Trimmomatic #
#----------------------------#
for i in $trim_path/*.err; do
    filename=$i
    # Get Sample Name
    grep "phred" $filename | cut -d " " -f3 | sed -e "s|$raw_file_path||g" | sed -e 's/\///g' | sed -e "s/$suffix//g" >> tmp_sample_name
    # Get Input Read Pairs
    grep "Input Read Pairs" $filename | cut -d " " -f4 >> tmp_input_read_pairs
    # Get Both Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f7 >> tmp_both_surviving
    # Get Both Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f8 | tr -d '()' >> tmp_both_surviving_percent
    # Get Forward Only Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f12 >> tmp_forward_surviving
    # Get Forward Only Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f13 | tr -d '()' >> tmp_forward_surviving_percent
    # Get Reverse Only Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f17 >> tmp_reverse_surviving
    # Get Reverse Only Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f18 | tr -d '()'>> tmp_reverse_surviving_percent
    # Get Dropped
    grep "Input Read Pairs" $filename | cut -d " " -f20 >> tmp_drop
    # Get Dropped Percent
    grep "Input Read Pairs" $filename | cut -d " " -f21 | tr -d '()' >> tmp_drop_percent;
done

paste -d "\t" tmp_sample_name tmp_input_read_pairs tmp_both_surviving tmp_both_surviving_percent tmp_forward_surviving tmp_forward_surviving_percent tmp_reverse_surviving tmp_reverse_surviving_percent tmp_drop tmp_drop_percent > tmp_trimmomatic_stats.tmp.txt
sort tmp_trimmomatic_stats.tmp.txt > trimmomatic_stats.table.txt

sed -i '1 i\Fastq_File\tInputReadPairs\tBothSurviving\tBothSurviving%\tForwardSurviving\tForwardSurviving%\tReverseSurviving\tReverseSurviving%\tReadPairsDropped\tReadPairsDropped%' trimmomatic_stats.table.txt
rm tmp_*

#------------------------------------------------------------------------------#

#---------------------#
# Get Stats From STAR #
#---------------------#

for i in $star_path/*.star.Log.final.out; do
    filename=$i
    # get the base fastq file name
    grep -l "Number of input reads" $filename | sed -e 's/.star.Log.final.out//g' | sed -e "s|$star_path||g" | sed -e 's/\///g' >> tmp_sample_name
    # Get Number of Input Reads
    grep "Number of input reads" $filename | cut -f2 >> tmp_input_reads
    # Get Uniquely Mapped Reads Number
    grep "Uniquely mapped reads number" $filename | cut -f2 >> tmp_unq_number
    # Get Uniquely Mapped Reads Percent
    grep "Uniquely mapped reads %" $filename | cut -f2 >> tmp_unq_percent
    # Get Multimapped Reads Number
    grep "Number of reads mapped to multiple loci" $filename | cut -f2 >> tmp_multi_number
    # Get Multimapped Reads Percent
    grep "% of reads mapped to multiple loci" $filename | cut -f2 >> tmp_multi_percent
    # Get Reads Mapped to Too Many Loci Number
    grep "Number of reads mapped to too many loci" $filename | cut -f2 >> tmp_too_many_loci_number
    # Get Reads Mapped to Too Many Loci Percent
    grep "% of reads mapped to too many loci" $filename | cut -f2 >> tmp_too_many_loci_percent
    # Get Percent Reads Unmapped: Too Many Mismatches
    grep "too many mismatches" $filename | cut -f2 >> tmp_too_many_mismatch_percent
    # Get Percent Reads Unmapped: Too Short
    grep "too short" $filename | cut -f2 >> tmp_too_short_percent
    # Get Percent Reads Unmapped: Other
    grep "% of reads unmapped: other" $filename | cut -f2 >> tmp_unmapped_other_percent;
done

paste -d "\t" tmp_sample_name tmp_input_reads tmp_unq_number tmp_unq_percent tmp_multi_number tmp_multi_percent tmp_too_many_loci_number tmp_too_many_loci_percent tmp_too_many_mismatch_percent tmp_too_short_percent tmp_unmapped_other_percent > tmp_star_stats.tmp.txt
sort tmp_star_stats.tmp.txt > star_stats.table.txt

sed -i '1 i\Fastq_File\tInputReadPairs\tUnqMap\tUnqMapPerc\tMultiMap\tMultiMapPerc\tTooManyLoci\tTooManyLociPerc\tUnmappedTooManyMismatch\tUnmappedTooShort\tUnmappedOther' star_stats.table.txt

rm tmp_*

#-----------------------------------------------------------------------------------------------#

#------------------------#
# Join the 3 stats files #
#------------------------#
# remove headers for merging
tail -n +2 fqc_stats.table.txt > fqc_stats.tmp.txt
tail -n +2 trimmomatic_stats.table.txt > trimmomatic_stats.tmp.txt
tail -n +2 star_stats.table.txt > star_stats.tmp.txt

file1=fqc_stats.tmp.txt
file2=trimmomatic_stats.tmp.txt
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10 <(sort -k1 $file1) <(sort -k1 $file2) > out1.txt

file1=out1.txt
file2=star_stats.tmp.txt
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11 <(sort -k1 $file1) <(sort -k1 $file2) > FQC_Trim_Aln_stats.tmp2.txt

rm *.tmp.txt
rm out1.txt

#-----------------------------------------------------------------------------------------------#

#-------------------------------------------------#
# Add Desired Sample name from metadata xlsx file #
#-------------------------------------------------#

# convert the xlsx file to txt tab separated file with R
./scripts/write_xlsx_to_txt.R -m $metadata

# join the above stats table with this new table (adding sample name)
file1=FQC_Trim_Aln_stats.tmp2.txt
file2=${metadata/xlsx/txt} 
join -j 1 -t $'\t' -o 1.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28 <(sort -k1 $file1) <(sort -k1 $file2) > FQC_Trim_Aln_stats.txt

# add heaaders to final table
sed -i '1 i\Fastq_File\tR1_raw_reads\tR2_raw_reads\tR1_read_length\tR2_read_length\tR1_GC%\tR2_GC%\tR1_dup%\tR2_dup%\tInputReadPairs\tBothSurviving\tBothSurviving%\tForwardSurviving\tForwardSurviving%\tReverseSurviving\tReverseSurviving%\tReadPairsDropped\tReadPairsDropped%\tInputReadPairs\tUnqMap\tUnqMapPerc\tMultiMap\tMultiMapPerc\tTooManyLoci\tTooManyLociPerc\tUnmappedTooManyMismatch\tUnmappedTooShort\tUnmappedOther' FQC_Trim_Aln_stats.txt

rm *tmp2.txt
