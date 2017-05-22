#!/bin/bash -l

mkdir COLO
mkdir COLO/20160128

cat COAD/20160128/LoF_Somatic.dat READ/20160128/LoF_Somatic.dat > COLO/20160128/LoF_Somatic.dat
cat COAD/20160128/Somatic_TCGA_VUS.txt READ/20160128/Somatic_TCGA_VUS.txt > COLO/20160128/Somatic_TCGA_VUS.txt
