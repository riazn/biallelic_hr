# Clear window
rm(list=ls())
library(data.table)

#     ______                 __  _
#    / ____/_  ______  _____/ /_(_)___  ____
#   / /_  / / / / __ \/ ___/ __/ / __ \/ __ \
#  / __/ / /_/ / / / / /__/ /_/ / /_/ / / / /
# /_/    \__,_/_/ /_/\___/\__/_/\____/_/ /_/

# Define Function
myMatrices_new <- function(Cancer,    #Cancer type TCGA abreviation
                       GeneList, #List of genes
                       LSTSCORE, #LST scores with sample Ids
                       Somatic,  #Somatic mutation file
                       Germline, #Germline mutation file 
                       LOH,      #Loh file 
                       CN,       #total CN file
                       SNP,      #If TRUE write down somMutsnps and germMutsnps file
                       score="LST")       #Default name of score is LST  
  {
    MDIR <- paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128",sep="");
    
    f <- GeneList 
    gene_list <-scan(file=f, what="character");
    nGenes <- length(gene_list);
    
    g <- LSTSCORE 
    LST_BRCA_Ids <- read.delim(g, sep="\t", header = F);
    Ids <- LST_BRCA_Ids$V2;
    nIds <- length(Ids);
    Ids <- as.character(Ids)
    Ids <-gsub(" ", "", Ids)
    
    N <- nGenes
    
    # Create binary matrix to store info Truncating somatic mutations
    somaticMut <- matrix(0, nIds, N) # initialize matrix with '0' in each entry; N == # of rows, nGenes == # of columns
    colnames(somaticMut) <- gene_list
    rownames(somaticMut) <- Ids
    
    # Read in mutation data (from sup table 6)
    f<-Somatic #paste(MDIR, Somatic, sep="");
    dat <- read.delim(file=f, header=F, stringsAsFactor=FALSE, sep = "\t")
    #Store info in matrix
    for (i in 1:dim(dat)[1]) {
      sample <- dat[i, "V1"];
      gene <- dat[i,"V2"];
      if(nchar(rownames(somaticMut)[1])==14){rownames(somaticMut)<-gsub(" ", "", rownames(somaticMut)); print("sucess_somatic")}
      if ( (gene %in% colnames(somaticMut)) && (sample %in% rownames(somaticMut)) ) { # If it finds the gene and the sample then enter '1' in somaticMut
        somaticMut[sample, gene] <- 1;
      } 
    }
    
    # Create matrix to store info LST
    LST <- matrix(0, nIds, 1) # initialize matrix with '0' in each entry; N == # of rows, nGenes == # of columns
    colnames(LST) <- "LST"
    rownames(LST) <- Ids
    
    #LST score
    f<-LSTSCORE 
    dat <- read.delim(file=f, header=F, stringsAsFactor=FALSE, sep = "\t")
    dat$V2 <-gsub(" ", "", Ids)
    #Store info in matrix
    for (i in 1:dim(dat)[1]) {
      lst <- dat[i, "V3"];
      sample <- gsub(" ","",as.character(dat[i, "V2"]));
      LST[sample,"LST"] <- lst;
    }
    
    # Create matrix to store info Germline Trunc Ruomu
    germlineMut <- matrix(0, nIds, N) # initialize matrix with '0' in each entry; N == # of rows, nGenes == # of columns
    colnames(germlineMut) <- gene_list
    rownames(germlineMut) <- Ids
    
    #Germline trunc muts
    f<-Germline 
    dat <- read.delim(file=f, header=F, stringsAsFactor=FALSE, sep = "\t")
    head(dat)
    dim(dat)
    
    #Store info in matrix
    for (i in 1:dim(dat)[1]) {
      sample<- dat[i, "V1"];
      gene <- dat[i, "V2"];
      if(nchar(rownames(germlineMut)[1])==14){rownames(germlineMut)<-gsub(" ", "", rownames(germlineMut)); print("sucess_germline")}
      if ( (gene %in% colnames(germlineMut)) && (sample %in% rownames(germlineMut)) ) { # If it finds the gene and the sample then enter '1' in germMut
        germlineMut[sample, gene] <- 1;
      } 
    }
    
    lohMut <- matrix(NA, nIds, N) # initialize matrix with '0' in each entry; N == # of rows, nGenes == # of columns
    colnames(lohMut) <- gene_list
    rownames(lohMut) <- Ids
    
    f<-LOH #paste(MDIR, LOH, sep="");
    dat <- read.delim(file=f, header=T, stringsAsFactor=FALSE,sep = "\t")
    rownames(dat) <- dat$TCGA_ID
    dat <- dat[,-1]
    
    #Store info in matrix: if the gene is absent from LOH matrix (no data available) introduce "NA"
    #For file dat of different length (less) than nIds
    if(nchar(rownames(dat)[1])==14){rownames(dat) <-gsub(" ", "", rownames(dat))}
    
    Ids_overlap <- rownames(LST)[which((Ids%in%rownames(dat)))]
    Ids_overlap <-gsub(" ", "", Ids_overlap)
    
    for(i in 1:length(gene_list)) {
      gene <- gene_list[i]
      if(gene%in%colnames(dat)){
        lohMut[Ids_overlap,gene] <- dat[Ids_overlap,gene]
      }else {
        lohMut[Ids_overlap,gene] <- rep("NA", length(Ids_overlap))  
      }
    }
    
    lohMut <-data.frame(lohMut)
    
    #TOTAL CN BRCA
    # Create matrix to store info TotalCNV
    cnvMut <- matrix(NA, nIds, N) # initialize matrix with '0' in each entry; N == # of rows, nGenes == # of columns
    colnames(cnvMut) <- gene_list
    rownames(cnvMut) <- Ids
    
    # Read values of TotalCNV in Matrix
    f<-CN 
    dat <- read.delim(file=f, header=T, stringsAsFactor=FALSE, sep = "\t")
    rownames(dat) <- dat$TCGA_ID
    dat <- dat[,-1]

    #For file dat of different length (less) than nIds
    if(nchar(rownames(dat)[1])==14){rownames(dat) <-gsub(" ", "", rownames(dat))}
    
    Ids_overlap <- rownames(LST)[which((Ids%in%rownames(dat)))]
    Ids_overlap <-gsub(" ", "", Ids_overlap)
    #Store info in matrix: if the gene is absent from CN matrix (no data available) introduce "NA"
    for(i in 1:length(gene_list)) {
      gene <- gene_list[i]
      if(gene%in%colnames(dat)){
        cnvMut[Ids_overlap,gene] <- dat[Ids_overlap,gene]
      }else {
        cnvMut[Ids_overlap,gene] <- rep("NA", length(Ids_overlap))  
      }
    }
    
    cnvMut <- data.frame(cnvMut)
    if(SNP){
      write.table(somaticMut, paste(MDIR, "/somaticMutsnps.csv", sep=""), sep="\t", quote = F, row.names = T, col.names = T);
      write.table(germlineMut, paste(MDIR, "/germlineMutsnps.csv", sep=""), sep="\t", quote = F, row.names = T, col.names = T);} else { 
      write.table(somaticMut, paste(MDIR, "/somaticMut.csv", sep=""), sep="\t", quote = F, row.names = T, col.names = T);
      write.table(germlineMut, paste(MDIR, "/germlineMut.csv", sep=""), sep="\t", quote = F, row.names = T, col.names = T)}

    write.table(LST, paste(MDIR, "/", sprintf("%s.csv",score), sep=""), sep="\t", row.names = T, quote=F, col.names = T)
    write.table(lohMut, paste(MDIR, "/lohMut.csv", sep=""), sep="\t", quote =F, row.names = T, col.names = T)
    write.table(cnvMut, paste(MDIR, "/cnvMut.csv", sep=""), sep="\t", quote = F,row.names = T, col.names = T)  
  }

#     __  ___      _
#    /  |/  /___ _(_)___
#   / /|_/ / __ `/ / __ \
#  / /  / / /_/ / / / / /
# /_/  /_/\__,_/_/_/ /_/

pan_cancer_types <- read.delim("Supplementary_Files/Cancers.txt",header=F)
pan_cancer_types <- gsub(" ","",as.character(pan_cancer_types$V1))
pan_cancer_types

##COAD and READ are merged under 'COLO'
pan_cancer_types <- pan_cancer_types[-c(5,24)]
pan_cancer_types <- c(pan_cancer_types,"COLO")
cancer <- "COLO"

##Read in Final_List of samples
samples_final <- read.delim("Supplementary_Files/Final_List.txt",header=F)
samples_final <- as.character(samples_final$V1)
length(samples_final)

## Re run matrix
Mast=list()

for (k  in 1:length(pan_cancer_types))
{
  
  Cancer <- as.character(pan_cancer_types[k])
  
  print(sprintf("%s-%i",Cancer,k))

  myMatrices_new(
    Cancer = Cancer,
    GeneList = "Supplementary_Files/Master_List.txt",
    LSTSCORE = sprintf("Supplementary_Files/LSTs/LST-%s.dat", Cancer),
    Somatic = paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/","LoF_Somatic.dat", sep=""), 
    Germline = sprintf("Supplementary_Files/Germline_files/LoF_germline_Ruomu_Git_%s_Id_Gene_R.txt",Cancer),
    LOH = sprintf("Supplementary_Files/LOH_calls/loh_%s_header.dat",Cancer),
    CN = sprintf("Supplementary_Files/CN_calls/totalcn_%s_header.dat",Cancer),
    SNP = F)

  myMatrices_new(
    Cancer = Cancer,
    GeneList = "Supplementary_Files/Master_List.txt",
    LSTSCORE = sprintf("Supplementary_Files/LSTs/LST-%s.dat", Cancer),
    Somatic =  paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/","Somatic_TCGA_VUS.txt", sep=""),
    Germline = sprintf("Supplementary_Files/Germline_files/LoF_germline_Ruomu_Git_%s_Id_Gene_R.txt",Cancer),
    LOH = sprintf("Supplementary_Files/LOH_calls/loh_%s_header.dat",Cancer),
    CN = sprintf("Supplementary_Files/CN_calls/totalcn_%s_header.dat",Cancer),
    SNP = T)
  
  somaticMutnonsyn <- read.table(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "somaticMutsnps.csv", sep=""), sep="\t", header = T)
  somaticMutnonsyn <- somaticMutnonsyn[which(as.character(rownames(somaticMutnonsyn))%in%samples_final),]

  somaticMut <- read.table(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "somaticMut.csv", sep=""), sep="\t", header = T)
  somaticMut <- somaticMut[which(as.character(rownames(somaticMut))%in%samples_final),]

  germlineMut <- read.table(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "germlineMut.csv", sep=""), sep="\t", header = T)
  germlineMut <- germlineMut[which(as.character(rownames(germlineMut))%in%samples_final),]

  LST <- read.table(paste("stddata__2016_01_28/",sprintf("%s/",Cancer), "20160128/", "LST.csv", sep=""), sep="\t", header = T)
  LST <- LST[which(as.character(rownames(LST))%in%samples_final),]

  lohMut <- read.table(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "lohMut.csv", sep=""), sep="\t", header = T)
  lohMut <- lohMut[which(as.character(rownames(lohMut))%in%samples_final),]

  cnvMut <- read.table(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "cnvMut.csv", sep=""), sep="\t", header = T)
  cnvMut <- cnvMut[which(as.character(rownames(cnvMut))%in%samples_final),]

  Somatic_Trunc <- read.delim(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/", "LoF_Somatic.dat", sep=""), header =F, sep ="\t")
  Somatic_Trunc <- Somatic_Trunc[which(as.character(Somatic_Trunc$V1)%in%samples_final),]

  Somatic_nonsyn <- read.delim(paste("stddata__2016_01_28/", sprintf("%s/",Cancer), "20160128/","Somatic_TCGA_VUS.txt", sep=""), header =F, sep ="\t")
  Somatic_nonsyn <- Somatic_nonsyn[which(as.character(Somatic_nonsyn$V1)%in%samples_final),]

  Mast[[k]] <- list(som_vus=somaticMutnonsyn,som_LoF=somaticMut,germ=germlineMut,lst=LST,loh=lohMut,cnv=cnvMut,Som_vus=Somatic_nonsyn,Som_LoF=Somatic_Trunc)
}

save(Mast,file = "Supplementary_Files/Matrices_all_Final_102genes.Rdat")
