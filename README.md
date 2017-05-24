# Pan-cancer analysis of bi-allelic alterations in homologous recombination (HR) DNA repair 
	The following provides step by step instructions to reproduce the key
	results from the manuscript above. Results can be reproduced
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
			- R version 3.2.1 
			- Python version 2.7.9 (please refer to Results_Figures_and_P_values/README.md for further details)
	This code has been tested on Mac OSX, Windows 7 and Linux environments
	Please email questions to riazn@mskcc.org

# A.) Download maf files from TCGA Firehose dated 2016/01/28: 
	- firehose_get (provided in the folder; make it an executable by typing 'chmod +x firehose_get'), 
	  which can be obtained from: 
	  https://confluence.broadinstitute.org/display/GDAC/Download

	- Use firehose_get command below to download maf files:

		`./firehose_get -tasks -b -only Mutation_Packager_Oncotated_Calls.Level_3 stddata 2016_01_28`

		Note: More cancer types data are downloaded than the 24 cancer types used in this paper.


# B.) Run the following scripts to prepare the data for analysis in R: 
	1. cd stddata__2016_01_28/
	2. Run Somatic_all_cancers_paper.sh
	3. Run merge_COAD_READ.sh to create the colo data
	4. Details of the code and output are provided within the README.md


# C.) Uncompress copy number data
	1. cd Supplementary_Files; tar -zxvf CopyNumberData.tgz

# D.) Run the following R codes in sequence to generate figures and P-values: 
  	1. Generate the mutation and allele specific copy number matrix:
		Prep_1_Mutation_Data_Generation.R 
		output --> Supplementary_Files/Matrices_all_Final_102genes.Rdat

	2. Generate the master matrix: 
		Prep_2_Matrix_Generation.R
		output --> Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt
		The different entries in the matrix are:
		0- wild type
		1- Germline biallelic pathogenic with LOH
		2- Germline biallelic pathogenic compound heterozygous (second somatic pathogenic hit)
		3- Germline biallelic VUS compound heterozygous (second somatic VUS hit)
		4- Germline mono allelic pahtogenic
		5- Somatic biallelic pathogenic with LOH
		6- Somatic biallelic pathogenic compund heterozygous (second somatic pathogenic hit)
		7- Somatic biallelic pathogenic VUS compund heterozygous (second somatic VUS hit)
		8- Somatic biallelic VUS with LOH
		9- Somatic monoallelic pathogenic
		10- Somatic monoallelic VUS

	3. Generate figures:
		- Incidence of biallelic pathogenic alterations across cancer types: 
		  Run Generate_Fig1a.R --> Results_Figures_and_P_Values/Fig1a.pdf
		- Incidence of biallelic pathogenic alterations per cancer type: 
		  Run Generate_Fig1b.R -->  >Results_Figures_and_P_Values/Fig1b.pdf
		- Association of LST and Signature3 in HBOC cancers:
		  Run Generate_Fig1c.R --> Results_Figures_and_P_Values/Fig1c ... .pdf
		- LST for each genotype:
		  Generate_Fig1d.R --> Results_Figures_and_P_Values/Fig1d.pdf
		- Association of LST and Signature3 across cancer types:
		  Generate_Fig1e.R --> Results_Figures_and_P_Values/Fig1e .... .pdf
		- Mutual Exclusivity analysis of biallelic alterations in HBOC cancers:
		  Run Generate_Fig2a_and_Input_for_WExT.R:
		  Generates --> Results_Figures_and_P_Values/Oncoprint_HBOCs_top16_LST15_OR_Dominant_MutSig3.txt
		  This result file will be fed to the cbio portal to generate the oncoprint:
				http://www.cbioportal.org/
				Use cbio portal --> tools --> Oncoprinter
		  Refer to Results_Figures_and_P_Values/README.md for further instructions on how to get the
		  mutual exclusivity test p  values using WExT (through the python codes provided: 
		  wext_run_comet_with_VUS.py and wext_run_comet_without_VUS.py).
 

	4. Generate supplementary figures:
		- Top 25 most frequently mutated HR-related genes:
		Generate_Suppl_Fig1a.R -- > Results_Figures_and_P_Values/Supplementary_Fig1a.pdf
		- Frequency of alterations per cancer type:
		Generate_Suppl_Fig1b.R -- > Results_Figures_and_P_Values/Supplementary_Fig1b.pdf
		- Breakdown of biallelic alteration types within the top 25 mutated genes:
		Generate_Suppl_Fig1c.R -- > Results_Figures_and_P_Values/Supplementary_Fig1c.pdf
		- Breakdown of biallelic alteration types across cancer types:
		Generate_Suppl_Fig1d.R -- > Results_Figures_and_P_Values/Supplementary_Fig1d.pdf

	Note: test_all_figures.R will run each of the figure scripts and save the results into the "Results_Figures_and_P_Values".
