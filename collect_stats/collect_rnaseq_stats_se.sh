#!/bin/bash

# Purpose: Collect stats from FastQC, Trimmomatic, STAR into one table
# For RNAseq data
# Currently for SE reads
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

##################
# User variables #
##################
# path to raw fastq files (which was used for trimmomatic)
# goal is to remove this string from the trimmomatic log, leaving only fastq file name
raw_file_path=/home/groups/hoolock2/u0/archive/Epigenetics_Core_rawdata/ECP67/220405_NS500681_0535_AHC7HCBGXK/RNA220304MS
# path to metadata excel file
metadata=/home/groups/hoolock2/u0/bd/Projects/ECP67/metadata/ECP67.metadata.xlsx
# Current directory is scripts directory
dir=$(pwd)
# change to project home dir
cd ..

#-----------#
# Variables #
#-----------#
# suffix of fastq files to remove
#suffix=_R1.fastq.gz
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

#------------------------------------------------------------------------

#-----------------------#
# Get Stats from FastQC #
#-----------------------#
cd $fqc_path
for f in `ls --color=none *.zip`; do unzip $f; done
# Grab Metrics
grep "Total Sequences" */fastqc_data.txt > tmp_reads.txt
grep "Sequence length" */fastqc_data.txt > tmp_readLength.txt
grep "Total Deduplicated" */fastqc_data.txt > tmp_duplicated.txt
grep "^%GC" */fastqc_data.txt > tmp_GC.txt
# Clean up Metrics
### transform the deduplicated percentage to duplicated percentage (by subtracting from 100)
awk -F'\t' '{print $1,(100-$2)}' tmp_duplicated.txt > tmp_duplicated2.txt
### Keep only the fastq prefix name in the first field of the tmp_reads.txt file
awk 'BEGIN {FS=OFS="\t"}{gsub("_R1.*","",$1)}1' tmp_reads.txt > tmp_file && mv tmp_file tmp_reads.txt
### Keep only the last field for remaining tmp files
cut -f2 tmp_readLength.txt > tmp2_readLength.txt
cut -f2 tmp_GC.txt > tmp2_GC.txt
cut -d ' ' -f4 tmp_duplicated2.txt > tmp2_duplicated.txt
# Combine FQC stats into one table
paste -d "\t" tmp_reads.txt tmp2_readLength.txt tmp2_GC.txt tmp2_duplicated.txt > tmp_fqc_stats.tmp.txt
sort tmp_fqc_stats.tmp.txt > fqc_stats.table.txt

# add a header
sed -i '1s/^/Fastq_File\traw_reads\tread_length\tGC%\tdup%\n/' fqc_stats.table.txt
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
    grep "phred" $filename | cut -d " " -f3 | sed -e "s|$raw_file_path||g" | sed -e 's/\///g' | awk 'BEGIN {FS=OFS="\t"}{gsub("_R1.*","",$1)}1' >> tmp_sample_name
    # Get Input Reads
    grep "Input Reads" $filename | cut -d " " -f3 >> tmp_input_read_pairs
    # Get Surviving
    grep "Input Reads" $filename | cut -d " " -f5 >> tmp_surviving
    # Get Surviving Percent
    grep "Input Reads" $filename | cut -d " " -f6 | tr -d '()' >> tmp_both_surviving_percent
    # Get Dropped
    grep "Input Reads" $filename | cut -d " " -f8 >> tmp_drop
    # Get Dropped Percent
    grep "Input Reads" $filename | cut -d " " -f9 | tr -d '()' >> tmp_drop_percent;
done

paste -d "\t" tmp_sample_name tmp_input_read_pairs tmp_surviving tmp_both_surviving_percent tmp_drop tmp_drop_percent > tmp_trimmomatic_stats.tmp.txt
sort tmp_trimmomatic_stats.tmp.txt > trimmomatic_stats.table.txt

sed -i '1 i\Fastq_File\tInputReads\tReadsSurviving\tSurviving%\tReadsDropped\tReadsDropped%' trimmomatic_stats.table.txt
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
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,2.2,2.3,2.4,2.5,2.6 <(sort -k1 $file1) <(sort -k1 $file2) > out1.txt

file1=out1.txt
file2=star_stats.tmp.txt
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11 <(sort -k1 $file1) <(sort -k1 $file2) > FQC_Trim_Aln_stats.tmp2.txt

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
join -j 1 -t $'\t' -o 1.1,2.2,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20 <(sort -k1 $file1) <(sort -k1 $file2) > FQC_Trim_Aln_stats.txt

# add heaaders to final table
sed -i '1 i\Fastq_File\tSample_Name\traw_reads\tread_length\tGC%\tdup%\tInputReads\tReeadsSurviving\tSurviving%\tReadsDropped\tReadsDropped%\tInputReads\tUnqMap\tUnqMapPerc\tMultiMap\tMultiMapPerc\tTooManyLoci\tTooManyLociPerc\tUnmappedTooManyMismatch\tUnmappedTooShort\tUnmappedOther' FQC_Trim_Aln_stats.txt

rm *tmp2.txt

# convert the full summary table txt file to an excel file
./scripts/write_txt_to_xlsx.R -i FQC_Trim_Aln_stats.txt
