# Create Somatic input files from maf files TCGA Firehose 28/01/16:  
    	1.) LoF_Somatic.dat (Loss of Function somatic mutations -LoF-) 
    	2.) Somatic_TCGA_VUS.txt (Somatic Variants of Unknown Significance -VUS-)
		^-- PEDRO:how are these created? (maybe this should be at end of README??)

# List of Cancer types
    Cancers.txt --> List of 24 cancer types


# Run main script from command line:
    bash Somatic_all_cancers_paper.sh (Windows, Linux)
    bash Somatic_all_cancers_paper_MACOS.sh ( Mac OSX)
    ###
    With this script we process the individual maf files per cancer type, merge them and 
    extract, all somatic calls that are LoF and those that are VUS.
    


	###Scripts contained into the main script from above:
	i.) Somatic_all_entries_paper.sh 
		merge, per cancer type, all individual sample mafs in a big single file
		
	ii.) Somatic_formatting_paper.sh 
		format the big single file from step above (per cancer type)
		
	iii.) Somatic_all_entries_VUS_LoF_paper.sh 
		Split the formated big single maf into the required input files: 
		LoF_Somatic.dat and Somatic_TCGA_VUS.txt


# Merge COAD and READ (COLO) in a single file 'colo-rectal':
  	merge_COAD_READ.sh
