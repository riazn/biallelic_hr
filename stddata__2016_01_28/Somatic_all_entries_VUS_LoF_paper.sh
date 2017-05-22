#!/bin/bash -l

awk '$10=="Splice_Site"' $1 > Splice.dat 
awk '$10=="Frame_Shift_Del"' $1 > Frame_Shift_Del.dat 
awk '$10=="Frame_Shift_Ins"' $1 > Frame_Shift_Ins.dat
awk '$10=="Nonstop_Mutation"' $1 > Nonstop.dat
awk '$10=="Nonsense_Mutation"' $1 > Nonsense.dat
awk '$10=="Start_Codon_Del"' $1 > Start_Codon_Del.dat
awk '$10=="Stop_Codon_Del"' $1 > Stop_Codon_Del.dat
awk '$10=="Start_Codon_Ins"' $1 > Start_Codon_Ins.dat
awk '$10=="Stop_Codon_Ins"' $1 > Stop_Codon_Ins.dat

## LoF somatic
cat Nonsense.dat Splice.dat Frame_Shift_Ins.dat Frame_Shift_Del.dat Nonstop.dat Start_Codon_Del.dat Stop_Codon_Del.dat Start_Codon_Ins.dat Stop_Codon_Ins.dat > LoF_Somatic.dat

rm Splice.dat
rm Frame_Shift_Ins.dat
rm Frame_Shift_Del.dat
rm Nonstop.dat
rm Nonsense.dat
rm Start_Codon_Del.dat
rm Stop_Codon_Del.dat
rm Start_Codon_Ins.dat
rm Stop_Codon_Ins.dat


## Variants of Unkown Significance somatic (VUS)
awk '$10 !="Silent" && $10 != "Intron" && $10 != "IGR" && $10 !~ /Flank/ && $10 != "RNA" && $10 !~ /UTR/ &&  $10 != "lincRNA"' $1 > Somatic_TCGA_nonsynonymous_LoF_VUS.txt
awk 'NR==FNR { a[$10]; next } !($10 in a)' LoF_Somatic.dat Somatic_TCGA_nonsynonymous_LoF_VUS.txt > Somatic_TCGA_VUS_transient.txt
sed '1,2d' Somatic_TCGA_VUS_transient.txt > Somatic_TCGA_VUS.txt
rm Somatic_TCGA_VUS_transient.txt
