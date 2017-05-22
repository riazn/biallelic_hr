#!/bin/bash -l

rm Somatic_TCGA_merged_all_entries_formatted_header.maf
out="Somatic_TCGA_merged_all_entries_formatted_header.maf"

awk '{if((NR >2 && $1!="Hugo_Symbol") || NR <= 2) print $0}' $1 | awk '{if($1 != "Unknown") print $0}' > $out
