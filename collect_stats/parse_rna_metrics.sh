

# Collect import percentages from output of CollectRnaSeqMetrics from picard tools

# grab file/sample name
grep --with-filename "PF_BASES" *.RNA_Metrics.txt | cut -d : -f1 > tmp_filename

# Grab column headers
#grep "PF_BASES" *.RNA_Metrics.txt | cut -f17,18,19,20,21 > tmp_headers

# Grab numbers
cat *.RNA_Metrics.txt | awk '/PF_BASES/ {getline;print}' | cut -f17,18,19,20,21 > tmp_numbers

#--------------------#

# combine results
paste -d "\t" tmp_filename tmp_numbers > summary_rna_metrics

# add column names
sed -i '1 i\FileName\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES' summary_rna_metrics

# remove tmp files
rm tmp_*
