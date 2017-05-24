#clear window: 
rm(list=ls())

### Figure 1D) Boxplot of bi-allelic pathogenic, bi-allelic vus, monoallelic pathogenic, monoallelic vus, rest for 
#HBOCs: Breast, Ovary and Prostate; Dominant Signature 3 cases.

getwd()
Score <- read.delim("Supplementary_Files/Scores_Combined_All_Cancers_allLST_allMutSigs_DominantSig.txt",header=T,sep=" ")
dim(Score)
mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat)

## Restrict to HBOCs

brca <- read.delim("Supplementary_Files/samples_looked_at//BRCA_Tumor_samples_Looked_at.csv",header=T,sep="\t")
ov <- read.delim("Supplementary_Files/samples_looked_at//OV_Tumor_samples_Looked_at.csv",header=T,sep="\t")
prad <- read.delim("Supplementary_Files/samples_looked_at//PRAD_Tumor_samples_Looked_at.csv",header=T,sep="\t")

hbocs <- rbind(brca,ov,prad)
dim(hbocs)

mat <- mat[,which(gsub("\\.","-",colnames(mat))%in%hbocs$ids)]
dim(mat)

##

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))

## 1,2,5,6 --> biallelic germline/somatic pathogenic

all_lst_1 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_biallelic_pathogenic <- colnames(mat)[which(mat[Gene,]==1 | mat[Gene,]==2 | mat[Gene,]==5 | mat[Gene,]==6)];
  length(mat_biallelic_pathogenic);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character((mat_biallelic_pathogenic))))]
  all_lst_1 <- c(all_lst_1,gsub(" ","",as.character(lst_samples)))
}


all_1 <- unique(all_lst_1)
length(all_1)
all_pathogenic <- all_1
##Annotate

samples_patho <- Score[which(as.character(Score$Id)%in%all_1),]
dim(samples_patho)
head(samples_patho)

samples_patho$Group <- "Biallelic Pathogenic"
head(samples_patho)
getwd()

## 3,7,8--> biallelic vus

all_lst_1 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_biallelic_pathogenic <- colnames(mat)[which(mat[Gene,]==3 | mat[Gene,]==7 | mat[Gene,]==8)];
  length(mat_biallelic_pathogenic);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character((mat_biallelic_pathogenic))))]
  all_lst_1 <- c(all_lst_1,gsub(" ","",as.character(lst_samples)))
}

length(all_lst_1)
all_1 <- unique(all_lst_1)
length(all_1)
all_1 <- all_1[which(!(all_1%in%as.character(samples_patho$Id)))]
all_biallelic_vus <- all_1
length(all_biallelic_vus)
##Annotate

samples_bi_vus <- Score[which(as.character(Score$Id)%in%all_biallelic_vus),]
head(samples_bi_vus)

samples_bi_vus$Group <- "Biallelic VUS"


## Mono allelic pathogenic mat --> 4,9
all_mono_patho <- NULL

for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_rest<- colnames(mat)[which(mat[Gene,]==4 | mat[Gene,]==9)];
  dim(mat_rest);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_rest)))]
  all_mono_patho <- c(all_mono_patho,gsub(" ","",as.character(lst_samples)))
}

length(all_mono_patho)

all_mono_patho <- unique(all_mono_patho)
all_mono_patho <- all_mono_patho[which(!(all_mono_patho%in%as.character(samples_patho$Id)))]
all_mono_patho <- all_mono_patho[which(!(all_mono_patho%in%as.character(samples_bi_vus$Id)))]
length(all_mono_patho)
samples_mono_patho <- Score[which(as.character(Score$Id)%in%all_mono_patho),]
head(samples_mono_patho)
dim(samples_mono_patho)
head(samples_mono_patho)
samples_mono_patho$Group <- "Mono allelic Pathogenic"

### Monoallelic VUS mat --> 10
all_mono_vus<- NULL

for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_rest<- colnames(mat)[which(mat[Gene,]==10)];
  dim(mat_rest);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_rest)))]
  all_mono_vus <- c(all_mono_vus,gsub(" ","",as.character(lst_samples)))
}

length(all_mono_vus)

all_mono_vus<- unique(all_mono_vus)
all_mono_vus <- all_mono_vus[which(!(all_mono_vus%in%as.character(samples_patho$Id)))]
all_mono_vus <- all_mono_vus[which(!(all_mono_vus%in%all_biallelic_vus))]
all_mono_vus <- all_mono_vus[which(!(all_mono_vus%in%all_mono_patho))]

samples_mono_vus <- Score[which(as.character(Score$Id)%in%all_mono_vus),]
dim(samples_mono_vus)

samples_mono_vus$Group <- "Mono allelic VUS"

### Rest mat --> 0
all_rest<- NULL

for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_rest<- colnames(mat)[which(mat[Gene,]==0)];
  dim(mat_rest);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_rest)))]
  all_rest <- c(all_rest,gsub(" ","",as.character(lst_samples)))
}

length(all_rest)

all_rest<- unique(all_rest)
all_rest <- all_rest[which(!(all_rest%in%as.character(samples_patho$Id)))]
all_rest <- all_rest[which(!(all_rest%in%all_biallelic_vus))]
all_rest<- all_rest[which(!(all_rest%in%all_mono_patho))]
all_rest<- all_rest[which(!(all_rest%in%all_mono_vus))]


samples_rest <- Score[which(as.character(Score$Id)%in%all_rest),]
head(samples_rest)
dim(samples_rest)
samples_rest$Group <- "Wild Type"

##Boxplot

nmc1 <- list(samples_patho$LST[which(samples_patho$dominant%in%"Signature.3")],samples_bi_vus$LST[which(samples_bi_vus$dominant%in%"Signature.3")],samples_mono_patho$LST[which(samples_mono_patho$dominant%in%"Signature.3")],samples_mono_vus$LST[which(samples_mono_vus$dominant%in%"Signature.3")],samples_rest$LST[which(samples_rest$dominant%in%"Signature.3")])
Boxplot1 <- sapply(nmc1,'[',seq(max(sapply(nmc1,length))))
colnames <- c("Biallelic_Patho","Biallelic_VUS","Monoallelic_Patho","Monoallelic_VUS","Rest")
colnames(Boxplot1) <- colnames 
Boxplot1 <- data.frame(Boxplot1)

nmc2 <- list(samples_patho$LST[which(!(samples_patho$dominant%in%"Signature.3"))],samples_bi_vus$LST[which(!(samples_bi_vus$dominant%in%"Signature.3"))],samples_mono_patho$LST[which(!(samples_mono_patho$dominant%in%"Signature.3"))],samples_mono_vus$LST[which(!(samples_mono_vus$dominant%in%"Signature.3"))],samples_rest$LST[which(!(samples_rest$dominant%in%"Signature.3"))])
Boxplot2 <- sapply(nmc2,'[',seq(max(sapply(nmc2,length))))
colnames <- c("Biallelic_Patho","Biallelic_VUS","Monoallelic_Patho","Monoallelic_VUS","Rest")
colnames(Boxplot2) <- colnames 
Boxplot2 <- data.frame(Boxplot2)

nmc3 <- list(samples_patho$LST,samples_bi_vus$LST,samples_mono_patho$LST,samples_mono_vus$LST,samples_rest$LST)
Boxplot3 <- sapply(nmc3,'[',seq(max(sapply(nmc3,length))))
colnames <- c("Biallelic_Patho","Biallelic_VUS","Monoallelic_Patho","Monoallelic_VUS","Rest")
colnames(Boxplot3) <- colnames 
Boxplot3 <- data.frame(Boxplot3)

pdf('Results_Figures_and_P_Values/Fig1d.pdf')
boxplot(Boxplot3 ,ylab = sprintf("%s","LST"), cex.sub=0.75, outline = F)
stripchart(Boxplot2,vertical = T, method = "jitter", pch = 1, col = 1, cex= 1, add = T)
stripchart(Boxplot1,vertical = T, method = "jitter", pch = 20, col = 2, cex= 1, add = T)
dev.off()

