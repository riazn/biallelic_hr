# Clear workspace
rm(list=ls())

library(data.table)

# Read in cancer types
pan_cancer_types <- read.delim("Supplementary_Files/Cancers.txt",header=F)
pan_cancer_types <- gsub(" ","",as.character(pan_cancer_types$V1))
pan_cancer_types

## COAD and READ are merged under 'COLO'
pan_cancer_types <- pan_cancer_types[-c(5,24)]
pan_cancer_types <- c(pan_cancer_types,"COLO")

## Read in curated list of samples

samples_final <- read.delim("Supplementary_Files/Final_List.txt",header=F)
samples_final <- as.character(samples_final$V1)
length(samples_final)

## Load TCGA mutation matrices
load("Supplementary_Files/Matrices_all_Final_102genes.Rdat")

## Get gene list
gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)


## Create the matrix with all the Biallelic pathogenic, Biallelic VUS, monoallelic pathogenic and monoallelic VUS calls for all cancer types.

## Initiate matrix

mat <- matrix(0,length(gene_list_hr),length(samples_final))
colnames(mat) <- samples_final
rownames(mat) <- gene_list_hr

dim(mat) 

List_som_snps <- NULL; 

for (i in 1:length(pan_cancer_types))
  {
    ##Loop over all cancers
    Cancer <- as.character(pan_cancer_types[i])
    
    for (k in 1:length(gene_list_hr))
      {
        ##Loop over all genes
         
        Gene <- gene_list_hr[k];
         
          
        print(sprintf("%s-%i",Cancer,k))
        
        
        somaticMutsnps <- Mast[[i]]$som_vus
        somaticMut <- Mast[[i]]$som_LoF
        germlineMut<- Mast[[i]]$germ
        lohMut <- Mast[[i]]$loh
        LST <- Mast[[i]]$lst
      
        Somatic_Trunc <- Mast[[i]]$Som_LoF
        Somaticsnps <- Mast[[i]]$Som_vus
        Somatic <- Somaticsnps[Somaticsnps$V2==Gene,]

        somatic_snps <- Somatic
        somatic_snps_Ids <- as.character(somatic_snps$V1)
        somatic_snps_Ids <- gsub(" ","",somatic_snps_Ids)
        
        somaticMutsnps_snps <- somaticMutsnps[which(rownames(somaticMutsnps)%in%somatic_snps_Ids),]
        lohMut_snps <- lohMut[which(rownames(lohMut)%in%somatic_snps_Ids),]
        
        #germline pathogenic with LOH: entry in mat --> 1
        
        germline <- rownames(lohMut)[which(germlineMut[,Gene]==1 & lohMut[,Gene]==1 & !is.na(lohMut[,Gene]))]
        germline <- germline[which(germline%in%colnames(mat))]
        mat[Gene,as.character(germline)] <- 1
        
        
        #germline biallelic pathogenic with a second somatic pathogenic hit: entry in mat --> 2
        
        germline_double <- rownames(lohMut)[which(germlineMut[,Gene]==1 & somaticMut[,Gene]==1)]
        germline_double <- germline_double[which(!(as.character(germline_double)%in%as.character(germline)))]
        germline_double <- germline_double[which(germline_double%in%colnames(mat))]
        germline_double_pathogenic <- germline_double
        mat[Gene,as.character(germline_double_pathogenic)] <- 2
        
        #germline biallelic with a second somatic VUS hit: entry mat --> 3
        
        germline_double <- rownames(lohMut)[which(germlineMut[,Gene]==1 & somaticMutsnps[,Gene]==1)]
        germline_double <- germline_double[which(!(as.character(germline_double)%in%as.character(germline)))]
        germline_double <- germline_double[which(!(as.character(germline_double)%in%as.character(germline_double_pathogenic)))]
        germline_double <- germline_double[which(germline_double%in%colnames(mat))]
        germline_double_VUS <- germline_double
        mat[Gene,as.character(germline_double_VUS)] <- 3

        
        #germline monoallelic pathogenic: entry in mat --> 4
        
        germline_mono <- rownames(lohMut)[which(germlineMut[,Gene]==1 & lohMut[,Gene]==0 & !is.na(lohMut[,Gene]))]
        germline_mono <- germline_mono[which(!(as.character(germline_mono)%in%as.character(germline_double_pathogenic)))]
        germline_mono <- germline_mono[which(!(as.character(germline_mono)%in%as.character(germline_double_VUS)))]
        germline_mono <- germline_mono[which(as.character(germline_mono)%in%as.character(colnames(mat)))]
        mat[Gene,as.character(germline_mono)] <- 4

        
        #Concatenate all germline hits
        germline_all <- c(as.character(germline),as.character(germline_double_pathogenic),as.character(germline_double_VUS),as.character(germline_mono))
        
        ##Somatic biallelic pathogenic with LOH mat <-- 5
        somatic_loh <- rownames(lohMut)[which(somaticMut[,Gene]==1 & lohMut[,Gene]==1 & !is.na(lohMut[,Gene]))]
        somatic_loh <- somatic_loh[which(as.character(somatic_loh)%in%as.character(colnames(mat)))]
        mat[Gene,as.character(somatic_loh)] <- 5
        
        ## Somatic biallelic pathogenic with a second somatic pathogenic hit: double mat <-- 6
       
        somatic_double <- rownames(lohMut)[which(somaticMut[,Gene]==1)]
        somatic_double <- somatic_double[which(!(as.character(somatic_double)%in%as.character(somatic_loh)))]
        somatic_double <- somatic_double[which(!(as.character(somatic_double)%in%as.character(germline_all)))]
     
        res <- Somatic_Trunc[which(as.character(Somatic_Trunc$V1)%in%somatic_double & as.character(Somatic_Trunc$V2)%in%Gene),]
        res <- data.table(res)
        setkey(res)
        res <- res[duplicated(res[list(V1, V2, V6),allow.cartesian=TRUE]) | duplicated(res[list(V1, V2, V6,V10),allow.cartesian=TRUE],fromLast = T)] #, as.character(V7), V8, V10, V12, V14, V43), allow.cartesian=TRUE])]
        res2 <- res[unique(res[list(V1,V2,V6,as.character(V7),V8,V10,V12,V14,V43),allow.cartesian=TRUE])]
        
        somatic_double_pathogenic <- NULL
        
        if(dim(res2)[1] !=0 & dim(res2)[1]==nrow(res)){
        res_path <- data.frame(res$V1,Cancer,Gene,res$V6,res$V7,res$V8,res$V10,res$V12,res$V14,res$V43)
        somatic_double_pathogenic <- res_path[,1]
        somatic_double_pathogenic <- somatic_double_pathogenic[which(as.character(somatic_double_pathogenic)%in%colnames(mat))]
        mat[Gene,as.character(somatic_double_pathogenic)] <- 6
        }
    
    ## Somatic biallelic with a second somatic VUS hit: mat <-- 7
    
    res <- Somatic[which(as.character(Somatic$V1)%in%somatic_double),]
    res <- data.table(res)
    setkey(res)
    res <- res[unique(res[list(V1, V2, V6, as.character(V7), V8, V10, V12, V14, V43), allow.cartesian=TRUE])]
    somatic_double <- res$V1
    res <- res[!duplicated(res$V1),]
    
    somatic_double_VUS <- NULL
    
    if(dim(res)[1]>0){
      
      res_vus <- data.frame(res$V1,Cancer,Gene,res$V3,res$V4,res$V5,res$V6,res$V7,res$V8,res$V9)
      somatic_double_VUS <- res_vus[,1]
      somatic_double_VUS <- somatic_double_VUS[which(as.character(somatic_double_VUS)%in%colnames(mat))] 
      somatic_double_VUS <- somatic_double_VUS[which(!(as.character(somatic_double_VUS)%in%as.character(somatic_double_pathogenic)))]
      mat[Gene,as.character(somatic_double_VUS)] <- 7
    }
   
    
    ### Somatic Biallelic VUS with LOH: mat <--8
    somaticsnps <- rownames(lohMut)[which(somaticMutsnps[,Gene]==1 & lohMut[,Gene]==1 & !is.na(lohMut[,Gene]))]
    res<- Somaticsnps[which(as.character(Somaticsnps$V1)%in%as.character(somaticsnps) & as.character(Somaticsnps$V2)%in%Gene),]
    res <- res[which(!(res$V1%in%germline_all)),]
    res <- res[which(res$V1%in%colnames(mat)),]
    dim(res)
    
    
    ## Annotate Biallelic VUS's
    if(dim(res)[1]>0){
      res <- data.frame(Id=res$V1,Cancer=Cancer,Gene=Gene,Chromosome=res$V6,Start=res$V7,End=res$V8,Type=res$V10,Ref_Allele=res$V12,Alt_Allele=res$V14,AAChange=res$V43)
      List_som_snps <-  rbind(List_som_snps,res)
    }
   
    if(dim(res)[1]>0){
    mat[Gene,as.character(res$Id)] <- 8;} 
   
    ## Somatic monoallelic pathogenic: mat --> 9
    somatic_mono <- rownames(lohMut)[which(somaticMut[,Gene]==1 & lohMut[,Gene]==0 & !is.na(lohMut[,Gene]))]
    somatic_mono <- somatic_mono[which(!(somatic_mono%in%germline_all))]
    somatic_mono <- somatic_mono[which(!(somatic_mono%in%somatic_double_pathogenic))]
    somatic_mono <- somatic_mono[which(!(somatic_mono%in%somatic_double_VUS))]
    somatic_mono <- somatic_mono[which(somatic_mono%in%colnames(mat))]
    mat[Gene,as.character(somatic_mono)] <- 9
    
    #somatic monoallelic VUS: mat --> 10
    somatic_mono_VUS <- rownames(somaticMutsnps_snps)[which(somaticMutsnps_snps[,Gene]==1 & lohMut_snps[,Gene]==0 & !is.null(lohMut_snps[,Gene]))]
    somatic_mono_VUS <- somatic_mono_VUS[which(!(somatic_mono_VUS%in%somatic_double_pathogenic))]
    somatic_mono_VUS <- somatic_mono_VUS[which(!(somatic_mono_VUS%in%germline_all))]
    somatic_mono_VUS <- somatic_mono_VUS[which(!(somatic_mono_VUS%in%somatic_double_VUS))] 
    somatic_mono_VUS <- somatic_mono_VUS[which(!(somatic_mono_VUS%in%somatic_mono))]
    somatic_mono_VUS <- somatic_mono_VUS[which(somatic_mono_VUS%in%colnames(mat))]
    
    mat[Gene,as.character(somatic_mono_VUS)] <- 10
    
    }
}

dim(mat)

write.table(mat,"Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",col.names=T,row.names=T,quote = F,sep="\t")

