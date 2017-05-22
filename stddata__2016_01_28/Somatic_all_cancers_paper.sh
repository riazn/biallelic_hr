#!/bin/bash -l

size=$(wc -l Cancers.txt | awk '{print $1}')

for k in `seq 1 $size`; do 

cancer=$(awk -v var=$k '{if(NR==var) print $1}' Cancers.txt)

cd $cancer/

echo $cancer

tar -zxvf gdac.broadinstitute.org_${cancer}.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0.tar.gz

bash ../Somatic_all_entries_paper.sh gdac.broadinstitute.org_${cancer}.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0 

bash ../Somatic_formatting_paper.sh Somatic_TCGA_merged_all_entries.maf 

bash ../Somatic_all_entries_VUS_LoF_paper.sh Somatic_TCGA_merged_all_entries_formatted_header.maf 

cd ../

done


