#clear window: 
rm(list=ls())

pan_cancer_types <- read.delim("Supplementary_Files/Cancers.txt",header=F)
pan_cancer_types <- gsub(" ","",as.character(pan_cancer_types$V1))
pan_cancer_types

##COAD and READ are merged under 'COLO'

pan_cancer_types <- pan_cancer_types[-c(5,24)]
pan_cancer_types <- c(pan_cancer_types,"COLO")
pan_cancer_types

## Upload genes
gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)

## load matrix

mat_cancers <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat_cancers)


## Patho bi and mono allelic hits per cancer accross genes: mat entries 1,2,5,6 (biallelic pathogenic) and 4,9 (monoallelic pathogenic)

List_All_Genes_HR=NULL
for (i in 1:length(pan_cancer_types))
{  
  Cancer <- as.character(pan_cancer_types[i]);
  print(Cancer);
  
  mat <- mat_cancers
  
  bi_patho <- NULL
  mono_patho <- NULL
  
  bams <- read.delim(sprintf("Supplementary_Files/samples_looked_at//%s_Tumor_samples_Looked_at.csv", Cancer),header=T,sep="\t")
  Samples_cancer <- bams$ids
  mat <- mat[,which(gsub("\\.","-",colnames(mat))%in%gsub(" ","",as.character(Samples_cancer)))]
  dim(mat)
  
  for(k in 1:length(gene_list_hr)){
    
    gene <- gene_list_hr[k];
    #print(gene);
    bi_ids_patho <- colnames(mat)[(which(mat[gene,]==1 | mat[gene,]==2 | mat[gene,]==5 | mat[gene,]==6))];
    mono_ids_patho <- colnames(mat)[(which(mat[gene,]==4 | mat[gene,]==9))];
    
    bi_patho <- c(bi_patho,bi_ids_patho);
    mono_patho <- c(mono_patho,mono_ids_patho);
  }
  bi_patho_pc <- length(unique(bi_patho))/length(Samples_cancer);
  mono_patho_pc <- length(unique(mono_patho))/length(Samples_cancer);
  
  if(is.null(List_All_Genes_HR)) {List_All_Genes_HR <- data.frame(biallelic_pathogenic=bi_patho_pc,monoallelic_pathogenic=mono_patho_pc)} else {
    List_All_Genes_HR <- rbind(List_All_Genes_HR,data.frame(biallelic_pathogenic=bi_patho_pc,monoallelic_pathogenic=mono_patho_pc))}
  
}


rownames(List_All_Genes_HR) <- pan_cancer_types
write.table(List_All_Genes_HR,"Results_Figures_and_P_values/Alterations_per_Cancer_across_genes_HR_Pathway_ordered_Biallelic_Monoallelic_Paper.txt", col.names=T,row.names=T,quote=F,sep="\t")



res <- List_All_Genes_HR
head(res)
res_count <- res[order(res$biallelic_pathogenic,decreasing = T),]
head(res_count)

table <- data.matrix(t(res_count))
dim(table)

rownames(table) <-  c("Biallelic Pathogenic","Monoallelic Pathogenic")


#Plot
counts <- table[,1:23]
outfile <- "Results_Figures_and_P_values/Supplementary_Fig1b.pdf"
pdf(outfile,width=10,height=7)
h <- barplot(counts,xaxt="n", main="",cex.main=0.9, ylab="Prevalence", cex.lab=1.5
             ,col=c("black","grey")
             ,beside=F)
labs <- as.vector(colnames(counts))
text(cex=1.5, x=h, y=-0.03,labs, xpd=TRUE, srt=90)
h<-legend("topright", legend=rownames(counts), title="Alteration Type:", fill=c("black","grey"),cex=1.5)
dev.off()


