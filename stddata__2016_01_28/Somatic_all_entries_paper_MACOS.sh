#!/bin/bash -l

# input is the directory containing all individual mafs for a give cancer type
# Example : Somatic_all_entries_paper.sh gdac.broadinstitute.org_BRCA.Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0 

rm Somatic_TCGA_merged_all_entries.maf

out="Somatic_TCGA_merged_all_entries.maf"

number=$(ls -l $1 | wc -l)

for k in $(seq 3 $number)
do

echo $k

file=$(ls -l $1 | awk -v VAR=$k '{if(NR==VAR) print $10}')

id=$(ls -l $1 | awk -v VAR=$k '{split($10,a,"."); split(a[1],b,"-"); if(NR==VAR) print b[1]"-"b[2]"-"b[3]}')

echo $id

cat $1/$file > dummy1.txt

#export ID=$id

sed '2,3d' dummy1.txt >> dummy2.txt

awk -v VAR=$id -F\t '{if(NR >1) print VAR"\t"$0}' dummy2.txt >> $out

rm dummy1.txt
rm dummy2.txt

done

