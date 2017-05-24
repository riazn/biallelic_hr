# This folder 'Results_Figures_and_P_Values' keeps all the final analysis results and figures.

# From Generate_Fig2a_and_Input_for_WExT.R script:
	1. It generates Fig2a.pdf
	2. It saves the Input for WExt in the Results_Figures_and_P_Values:
	   i) vec_top_16_genes_without_VUS.txt
	   ii) vec_top_16_genes_with_VUS.txt
	   iii) genes_top_16_genes_without_VUS.txt
	   iV) genes_top_16_genes_with_VUS.txt
	  
	3. Run the python scripts below to generate the p value of the mutual exclusivity test
		wext_run_comet_with_VUS.py
		wext_run_comet_without_VUS.py

		 - WExT source code and instructions for installation
		 		https://github.com/raphael-group/wext

		 - Mutual exclusivity test Without VUS
			time python wext_run_comet_without_VUS.py

		- Mutual exclusivity test With VUS
			time python wext_run_comet_with_VUS.py  
