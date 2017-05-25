# Pan-cancer analysis of bi-allelic alterations in homologous recombination (HR) DNA repair 
	The following provides step by step instructions to reproduce the key
	results from the manuscript above. Results can be reproduced
	from the very beginning by downloading MAF files from TCGA (start at
	step A below) or from a precomputed matrix of mutation and copy
	numbers calls by skipping to step (D.3)
	
	Pre-computed results included are: Large-scale transitions (LST) values
	have been provided as these require access to Affy SNP arrays (TCGA
	Level I) data access. Similarly the proportions of mutations due to
	signature 3 have already also been provided. Loss-of-heterozygosity and
	total copy number results are also provided as these require access to
	TCGA Level I data (SNP6 Arrays). 
	
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

		./firehose_get -tasks -b -only Mutation_Packager_Oncotated_Calls.Level_3 stddata 2016_01_28

		Note: More cancer types data are downloaded than the 24 cancer types used in this paper.
		
		


# B.) Run the following scripts to prepare the data for analysis in R: 
	1. cd stddata__2016_01_28/
	2. Run Somatic_all_cancers_paper.sh (Windows, Linux)
	   Run Somatic_all_cancers_paper_MACOS.sh (Mac OSX)
	   Note: ignore the 'rm Somatic_TCGA_merged_all_entries.maf: No such file or directory' warning
	   
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
		This creates an integer matrix where the columns represent TCGA cases, and the 
		rows, individual genes. The values correspond to:
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


	3. Generate figures Regarding incidence of bi-allelic and mono-allelic mutations:
		- Incidence of biallelic pathogenic alterations across cancer types: 
		  Run Generate_Fig1a.R --> Results_Figures_and_P_Values/Fig1a.pdf
		  
		- Incidence of biallelic pathogenic alterations per cancer type: 
		  Run Generate_Fig1b.R -->  >Results_Figures_and_P_Values/Fig1b.pdf
		  
		- Top 25 most frequently mutated HR-related genes:
		Generate_Suppl_Fig1a.R -- > Results_Figures_and_P_Values/Supplementary_Fig1a.pdf
		
		- Frequency of alterations per cancer type:
		Generate_Suppl_Fig1b.R -- > Results_Figures_and_P_Values/Supplementary_Fig1b.pdf
		
		- Breakdown of biallelic alteration types within the top 25 mutated genes:
		Generate_Suppl_Fig1c.R -- > Results_Figures_and_P_Values/Supplementary_Fig1c.pdf
		
		- Breakdown of biallelic alteration types across cancer types:
		Generate_Suppl_Fig1d.R -- > Results_Figures_and_P_Values/Supplementary_Fig1d.pdf

		  
	4. Evaluate the association between genomic evidence of HR deficiency
	   (LST & Mutational Signature 3) and bi-allelic genetic alterations in HR Genes
	  
		- Association of LST and Signature3 in HBOC cancers:
		  Run Generate_Fig1c.R --> Results_Figures_and_P_Values/Fig1c ... .pdf
		  Generate_Fig1d.R --> Results_Figures_and_P_Values/Fig1d.pdf
		  
		- Association of LST and Signature3 across cancer types:
		  Generate_Fig1e.R --> Results_Figures_and_P_Values/Fig1e .... .pdf

 
	5. Determine HR genes are mutually exclusive in HBOC cancers
		 - Run Generate_Fig2a_and_Input_for_WExT.R:
		 Generates --> Results_Figures_and_P_Values/Oncoprint_HBOCs_top16_LST15_OR_Dominant_MutSig3.txt

		 - To get the mutual exclusivity test p values first install WExT: 
	 	 https://github.com/raphael-group/wext

	         - Run Mutual exclusivity test Without VUS
		 time python wext_run_comet_without_VUS.py

		 - Mutual exclusivity test With VUS
		 time python wext_run_comet_with_VUS.py
		
		 ### Python codes provided here for the above 
		  		- wext_run_comet_with_VUS.py 
		  		- wext_run_comet_without_VUS.py

		  - An oncoprint can be generated from cbio portal:
				http://www.cbioportal.org/
				Use cbio portal --> tools --> Oncoprinter

	Note: test_all_figures.R will run each of the figure scripts and save the results into the "Results_Figures_and_P_Values".
