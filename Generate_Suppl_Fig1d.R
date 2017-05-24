rm(list=ls())

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)

pan_cancer_types <- read.delim("Supplementary_Files/Cancers.txt",header=F)
pan_cancer_types <- gsub(" ","",as.character(pan_cancer_types$V1))
pan_cancer_types

##COAD and READ are merged under 'COLO'

pan_cancer_types <- pan_cancer_types[-c(5,24)]
pan_cancer_types <- c(pan_cancer_types,"COLO")
pan_cancer_types
## load matrix

mat_cancers <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat_cancers)


## Per cancer accross genes mut types
List_All_Genes_HR=NULL
for (i in 1:length(pan_cancer_types))
{  
  Cancer <- as.character(pan_cancer_types[i]);
  print(Cancer);
  
  germ_bi_patho <- NULL
  germ_bi_patho_double <- NULL
  germ_bi_vus <- NULL
  
  som_bi_patho <- NULL
  som_bi_patho_double <- NULL
  som_bi_vus <- NULL
  
  som_loh <- NULL
  
  
  bams <- read.delim(sprintf("Supplementary_Files/samples_looked_at//%s_Tumor_samples_Looked_at.csv", Cancer),header=T,sep="\t")
  Samples_cancer <- bams$ids
  mat <- mat_cancers[,which(gsub("\\.","-",colnames(mat_cancers))%in%gsub(" ","",as.character(Samples_cancer)))]
  

for(k in 1:length(gene_list_hr)){
  
  
  gene <- gene_list_hr[k];

  germ_bi_patho_ids <- colnames(mat)[which(mat[gene,]==1)];
  germ_bi_patho_double_ids <- colnames(mat)[which(mat[gene,]==2)];
  germ_bi_vus_ids <- colnames(mat)[which(mat[gene,]==3)];
  
  som_bi_patho_ids <- colnames(mat)[which(mat[gene,]==5)];
  som_bi_patho_double_ids <- colnames(mat)[which(mat[gene,]==6)];
  som_bi_vus_ids <- colnames(mat)[which(mat[gene,]==7)];
  
  som_loh_ids <- colnames(mat)[which(mat[gene,]==8)];
  
  germ_bi_patho <- c(germ_bi_patho,germ_bi_patho_ids);
  germ_bi_patho_double <- c(germ_bi_patho_double,germ_bi_patho_double_ids);
  germ_bi_vus <- c(germ_bi_vus,germ_bi_vus_ids);
  
  som_bi_patho <- c(som_bi_patho, som_bi_patho_ids);
  som_bi_patho_double <- c(som_bi_patho_double, som_bi_patho_double_ids);
  som_bi_vus <- c(som_bi_vus, som_bi_vus_ids);
  
  som_loh <- c(som_loh, som_loh_ids);
  }

  res <- data.frame(germ_biallelic_pathogenic=length(unique(germ_bi_patho))/length(Samples_cancer),germ_biallelic_pathogenic_double=length(unique(germ_bi_patho_double))/length(Samples_cancer),germ_biallelic_vus=length(unique(germ_bi_vus))/length(Samples_cancer),som_biallelic_pathogenic=length(unique(som_bi_patho))/length(Samples_cancer),som_biallelic_pathogenic_double=length(unique(som_bi_patho_double))/length(Samples_cancer),som_biallelic_vus=length(unique(som_bi_vus))/length(Samples_cancer),som_loh=length(unique(som_loh))/length(Samples_cancer))
  if(is.null(List_All_Genes_HR)) {List_All_Genes_HR <- res} else {
    List_All_Genes_HR <- rbind(List_All_Genes_HR,res)}
}

rownames(List_All_Genes_HR) <- pan_cancer_types
write.table(List_All_Genes_HR,"Supplementary_Files/Alterations_per_Cancer_across_genes_HR_Pathway_ordered_Granular_Paper.txt", col.names=T,row.names=T,quote=F,sep="\t")



res <- List_All_Genes_HR
colnames(res) <- c("Germline Trunc LOH","Germline Trunc Somatic Trunc","Germline Trunc Somatic not Trunc","Somatic Trunc LOH","Somatic Trunc Somatic Trunc","Somatic Trunc Somatic not Trunc","Somatic not Trunc LOH")
head(res)

##Count pathogenic events
Count=NULL
for(k in 1:dim(res)[1])
{
  cancer <- rownames(res)[k]
  print(cancer)
  count <- res[cancer,"Germline Trunc LOH"]+res[cancer,"Germline Trunc Somatic Trunc"]+res[cancer,"Somatic Trunc LOH"]+res[cancer,"Somatic Trunc Somatic Trunc"]
  Count <- c(Count,count)
}

count <- data.frame(Count)
rownames(count) <- rownames(res)
res_count <- cbind(res,count)
res_count <- data.frame(res_count)
res_count <- res_count[order(res_count$Count,decreasing = T),]
head(res_count)

##Merge compund het counts, germ and som

res_merged <- data.frame(res_count$Germline.Trunc.LOH,res_count$Germline.Trunc.Somatic.Trunc + res_count$Germline.Trunc.Somatic.not.Trunc,res_count$Somatic.Trunc.LOH,res_count$Somatic.Trunc.Somatic.Trunc + res_count$Somatic.Trunc.Somatic.not.Trunc,res_count$Somatic.not.Trunc.LOH, res_count$Count)
head(res_merged)
rownames(res_merged) <- rownames(res_count)

table <- data.matrix(t(res_merged))
dim(table)
head(table)
rownames(table) <-  c("Germ Patho LOH","Germ Comp Het","Som Patho LOH","Som Comp Het","Som VUS LOH","Count")


#Plot
counts <- table[-6,]
outfile <- "Results_Figures_and_P_Values/Supplementary_Fig1d.pdf"
pdf(outfile,width=10,height=7)
h <- barplot(counts,xaxt="n", main="",cex.main=0.9, ylab="Prevalence", cex.lab=1.5
             ,col=c("black","blue","grey","lightblue","darkgreen")
             ,beside=F)
labs <- as.vector(colnames(counts))
text(cex=1.5, x=h, y=-0.025,labs, xpd=TRUE, srt=90)
h<-legend("topright", legend=rownames(counts), title="Alteration Type:", fill=c("black","blue","grey","lightblue","darkgreen"),cex=1)
dev.off()

