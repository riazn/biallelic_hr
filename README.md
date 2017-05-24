# Pan-cancer analysis of bi-allelic alterations in homologous recombination (HR) DNA repair 
	The following provides step by step instructions to reproduce key
	the results from the manuscript above. Results can be reproduced
	from the very beginning by downloading MAF files from TCGA (start at
	step A below) or from a precomputed matrix of mutation and copy
	numbers calls by skipping to step (D.)
	
	Pre-computed results included are: Large-scale transitions (LST) values
	have been provided as these require access to Affy SNP arrays (TCGA
	Level I) data access. Similarly the proportions of mutations due to
	signature 3 have already also been provided. Loss-of-heterozygosity and
	total copy number results are also provided as these require access to
	TCGA Level I data (SNP6 Arrays)
	
	Please note the following software is required:
			- R version 3.x.x
			- 
	This code has been tested on Mac OSX and Linux environments
	Please email questions to riazn@mskcc.org

# A.) Download maf files from TCGA Firehose dated 2016/01/28: 
	- firehose_get (provided in the folder), which can be obtained from: 
	  https://confluence.broadinstitute.org/display/GDAC/Download

	- Use firehose_get command below to download maf files:

		`./firehose_get -tasks -b -only Mutation_Packager_Oncotated_Calls.Level_3 stddata 2016_01_28`

		Note: More cancer types data are downloaded than the 24 cancer types used in this paper.


# B.) Run the following scripts to prepare the data for analysis in R: 
	1. Run stddata__2016_01_28/Somatic_all_cancers_paper.sh
	2. Run stddata__2016_01_28/merge_COAD_READ.sh to create the colo data


# C.) Uncompress copy number data
	1. cd Supplementary_Files; tar -zxvf CopyNumberData.tgz

# D.) Run the following R codes in sequence to generate figures and P-values: 

	1. Generate the mutation and allele specific copy number matrix:
		Prep_1_Mutation_Data_Generation.R 

	2. Generate the master matrix: 
		Prep_2_Matrix_Generation.R

	3. Generate figures:
		Generate_Fig1a.R --> Results_Figures_and_P_Values/Fig1a.pdf
		Generate_Fig1b.R --> Results_Figures_and_P_Values/Fig1b.pdf
		Generate_Fig1c.R --> Results_Figures_and_P_Values/Fig1c ... .pdf
		Generate_Fig1d.R --> Results_Figures_and_P_Values/Fig1d.pdf
		Generate_Fig1e.R --> Results_Figures_and_P_Values/Fig1e .... .pdf

		Generate_Fig2a.R --> Results_Figures_and_P_Values/Oncoprint_HBOCs_top16_LST15_OR_Dominant_MutSig3.txt
			This result file will be fed to the cbio portal to generate the oncoprint:
				http://www.cbioportal.org/
				Use cbio portal --> tools --> Oncoprinter
 

	4. Generate supplementary figures:
		Generate_Suppl_Fig1a.R -- > Results_Figures_and_P_Values/Supplementary_Fig1a.pdf
		Generate_Suppl_Fig1b.R -- > Results_Figures_and_P_Values/Supplementary_Fig1b.pdf
		Generate_Suppl_Fig1c.R -- > Results_Figures_and_P_Values/Supplementary_Fig1c.pdf
		Generate_Suppl_Fig1d.R -- > Results_Figures_and_P_Values/Supplementary_Fig1d.pdf

	Note: test_all_figures.R will run each of the figure scripts and save the results into the "Results_Figures_and_P_Values".

