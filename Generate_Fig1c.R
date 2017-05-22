#clear window: 
rm(list=ls())

### Figure 1C) Boxplot of bi-allelic pathogenic, bi-allelic vus, monoallelic pathogenic, monoallelic vus, rest for 
#HBOCs: Breast, Ovary and Prostate

Score <- read.delim("Supplementary_Files/Scores_Combined_BRCA_OV_PRAD_LSTall_MutSig3all.txt",header=T,sep="\t")
head(Score)
mat <- read.delim("Supplementary_Files/matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat)

mat_cancers <- mat[,which(gsub("\\.","-",colnames(mat))%in%gsub(" ","",as.character(Score$Id)))]
dim(mat_cancers)


gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))

mat <- mat_cancers

## biallelic pathogenic: entries 1,2,5 and 6

all_lst_1 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_biallelic_pathogenic <- colnames(mat)[which(mat[Gene,]==1 | mat[Gene,]==2 | mat[Gene,]==5 | mat[Gene,]==6)];
  dim(mat_biallelic_pathogenic);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character((mat_biallelic_pathogenic))))]
  all_lst_1 <- c(all_lst_1,gsub(" ","",as.character(lst_samples)))
}


all_lst_1 <- unique(all_lst_1)
lst_biallelic_pathogenic <- Score$LST[which(gsub(" ","",as.character(Score$Id))%in%unique(all_lst_1))]
length(lst_biallelic_pathogenic)
length(all_lst_1)

## biallelic VUS: entries 3,7,8

all_lst_2 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_biallelic_VUS <- colnames(mat)[which(mat[Gene,]==3 | mat[Gene,]==7 | mat[Gene,]==8)];
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character((mat_biallelic_VUS))))]
  all_lst_2 <- c(all_lst_2,gsub(" ","",as.character(lst_samples)))
}

length(all_lst_2)
all_lst_2 <- unique(all_lst_2)
uniques <- unique(all_lst_2[which(!(all_lst_2%in%all_lst_1))])
all_lst_2 <- uniques
lst_biallelic_vus<- Score$LST[which(gsub(" ","",as.character(Score$Id))%in%all_lst_2)]
length(lst_biallelic_vus)
length(all_lst_2)


## monoallelic pathogenic: entries 4,9

all_lst_3 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_mono_pathogenic <- colnames(mat)[which(mat[Gene,]==4 | mat[Gene,]==9)];
  dim(mat_mono_pathogenic);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_mono_pathogenic)))]
  all_lst_3 <- c(all_lst_3,gsub(" ","",as.character(lst_samples)))
}

uniques <- c(all_lst_1,all_lst_2)
uniques <- unique(all_lst_3[which(!(all_lst_3%in%uniques))])
all_lst_3 <- uniques
lst_monoallelic_pathogenic <- Score$LST[which(gsub(" ","",as.character(Score$Id))%in%all_lst_3)]
length(lst_monoallelic_pathogenic)
length(all_lst_3)

## monoallelic VUS: entry 10

all_lst_4 <- NULL
for(k in 1:length(gene_list_hr)){
  
  lst <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_mono_VUS <- colnames(mat)[which(mat[Gene,]==10)];
  dim(mat_mono_VUS);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_mono_VUS)))]
  all_lst_4 <- c(all_lst_4,gsub(" ","",as.character(lst_samples)))
}


uniques <- c(all_lst_1,all_lst_2,all_lst_3)
uniques <- unique(all_lst_4[which(!(all_lst_4%in%uniques))])
all_lst_4 <- uniques
length(all_lst_4)
lst_monoallelic_VUS <- Score$LST[which(gsub(" ","",as.character(Score$Id))%in%all_lst_4)]
length(lst_monoallelic_VUS)

##unaltered cases (res): entry 0

all_lst_rest <- NULL

for(k in 1:length(gene_list_hr)){
  
  lst_samples <- NULL
  Gene <- gene_list_hr[k];
  print(Gene);
  mat_rest<- colnames(mat)[which(mat[Gene,]==0)];
  dim(mat_rest);
  lst_samples <- Score$Id[which(gsub(" ","",as.character(Score$Id))%in%gsub("\\.","-",as.character(mat_rest)))]
  all_lst_rest <- c(all_lst_rest,gsub(" ","",as.character(lst_samples)))
}


list_all_samples <- c(all_lst_1,all_lst_2,all_lst_3,all_lst_4)
list_all_samples <- unique(list_all_samples)
length(list_all_samples)
lst_rest_samples <- all_lst_rest[which(!(all_lst_rest%in%list_all_samples))]
unique <- unique(lst_rest_samples)
length(unique)
all_rest <- unique
lst_rest <- Score$LST[which(gsub(" ","",as.character(Score$Id))%in%all_rest)]
length(lst_rest)
length(unique)

##BoxPlot LST

res <- c(lst_biallelic_pathogenic,lst_biallelic_vus,lst_monoallelic_pathogenic,lst_monoallelic_VUS,lst_rest)
length(res)

nmc <- list(lst_biallelic_pathogenic,lst_biallelic_vus,lst_monoallelic_pathogenic,lst_monoallelic_VUS,lst_rest)
Boxplot <- sapply(nmc,'[',seq(max(sapply(nmc,length))))
colnames <- c("Biallelic_Patho","Biallelic_VUS","Monoallelic_Patho","Monoallelic_VUS","Rest")
colnames(Boxplot) <- colnames 
Boxplot <- data.frame(Boxplot)

Pval <- NULL
pval <- wilcox.test(Boxplot$Biallelic_Patho,Boxplot$Rest)$p.val
Pval <-c(Pval,pval)
pval <- wilcox.test(Boxplot$Biallelic_VUS,Boxplot$Rest)$p.val
Pval <-c(Pval,pval)
pval <- wilcox.test(Boxplot$Monoallelic_Patho,Boxplot$Rest)$p.val
Pval <-c(Pval,pval)
pval <- wilcox.test(Boxplot$Monoallelic_VUS,Boxplot$Rest)$p.val
Pval <-c(Pval,pval)

##Output Pval
print(Pval)
write.table(Pval,"Results_Figures_and_P_values/P_values_LST_Fig1c.txt",col.names=F,row.names=F,quote=F,sep="\t")

pdf('Results_Figures_and_P_values/Fig1c_LST.pdf')
boxplot(Boxplot ,ylab = sprintf("%s","LST"), cex.sub=0.75, outline = F)
stripchart(Boxplot,vertical = T, method = "jitter", pch = 1, col = 1, cex= 1, add = T)
dev.off()

### Bosplot Signature 3

## 1 --> biallelic pathogenic

m3_biallelic_pathogenic <- Score$M3[which(gsub(" ","",as.character(Score$Id))%in%all_lst_1)]
length(m3_biallelic_pathogenic)

## 2 --> biallelic VUS

m3_biallelic_vus <- Score$M3[which(gsub(" ","",as.character(Score$Id))%in%all_lst_2)]
length(m3_biallelic_vus)


## 3 --> monoallelic pathogenic

m3_monoallelic_pathogenic <- Score$M3[which(gsub(" ","",as.character(Score$Id))%in%all_lst_3)]
length(m3_monoallelic_pathogenic)



## 4 --> monoallelic VUS

m3_monoallelic_vus <- Score$M3[which(gsub(" ","",as.character(Score$Id))%in%all_lst_4)]
length(m3_monoallelic_vus)

## 0 --> unaltered cases

m3_rest <- Score$M3[which(gsub(" ","",as.character(Score$Id))%in%unique(all_rest))]
length(m3_rest)

##BoxPlot

nmc_m3 <- list(m3_biallelic_pathogenic,m3_biallelic_vus,m3_monoallelic_pathogenic,m3_monoallelic_vus,m3_rest)
Boxplot_m3 <- sapply(nmc_m3,'[',seq(max(sapply(nmc_m3,length))))
colnames <- c("Biallelic_Patho","Biallelic_VUS","Monoallelic_Patho","Monoallelic_VUS","Rest")
colnames(Boxplot_m3) <- colnames 
Boxplot_m3 <- data.frame(Boxplot_m3)


Pval <- NULL
pval_lst <- wilcox.test(Boxplot_m3$Biallelic_Patho,Boxplot_m3$Rest,alternative = "greater")$p.val
Pval <-c(Pval,pval_lst)
pval_lst <- wilcox.test(Boxplot_m3$Biallelic_VUS,Boxplot_m3$Rest,alternative = "greater")$p.val
Pval <-c(Pval,pval_lst)
pval_lst <- wilcox.test(Boxplot_m3$Monoallelic_Patho,Boxplot_m3$Rest,alternative = "greater")$p.val
Pval <-c(Pval,pval_lst)
pval_lst <- wilcox.test(Boxplot_m3$Monoallelic_VUS,Boxplot_m3$Rest,alternative = "greater")$p.val
Pval <-c(Pval,pval_lst)

print(Pval)
##Output Pval

write.table(Pval,"Results_Figures_and_P_values/P_values_Signature_3_Fig1c.txt",col.names=F,row.names=F,quote=F,sep="\t")

pdf('Results_Figures_and_P_values/Fig1c_Signature3.pdf')
boxplot(Boxplot_m3 ,ylab = sprintf("%s","Signature 3"), cex.sub=0.75, outline = F)
stripchart(Boxplot_m3,vertical = T, method = "jitter", pch = 1, col = 1, cex= 1, add = T)
dev.off()

