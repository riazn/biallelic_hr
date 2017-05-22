#clear window: 
rm(list=ls())

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)

## load matrix

mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat)


## Patho bi and mono allelic hits per gene accross cancers: mat entries 1,2,3,5,6,7,8 (biallelic pathogenic) and 4,9 (monoallelic pathogenic)

List_All_Genes_HR=NULL

for(k in 1:length(gene_list_hr)){
  
  bi_patho <- 0
  mono_patho <- 0

  
  gene <- gene_list_hr[k];
  print(gene);
  bi_num_patho <- length(which(mat[gene,]==1 | mat[gene,]==2 | mat[gene,]==5 | mat[gene,]==6));
  mono_num_patho <- length(which(mat[gene,]==4 | mat[gene,]==9));
  
  bi_patho <- bi_patho + bi_num_patho;
  mono_patho <- mono_patho + mono_num_patho;

  
  if(is.null(List_All_Genes_HR)) {List_All_Genes_HR <- data.frame(biallelic_pathogenic=bi_patho,monoallelic_pathogenic=mono_patho)} else {
   List_All_Genes_HR <- rbind(List_All_Genes_HR,data.frame(biallelic_pathogenic=bi_patho,monoallelic_pathogenic=mono_patho))}

}


rownames(List_All_Genes_HR) <- gene_list_hr

res <- List_All_Genes_HR
head(res)
res_count <- res[order(res$biallelic_pathogenic,decreasing = T),]
head(res_count)

table <- data.matrix(t(res_count))
dim(table)

rownames(table) <-  c("Biallelic Pathogenic","Monoallelic Pathogenic")


#Plot
counts <- table[,1:25]

pdf('Results_Figures_and_P_values/Supplementary_Fig1a.pdf')
h <- barplot(counts,xaxt="n", main="",cex.main=0.9, ylab="Altered Tumors", cex.lab=1.5
             ,col=c("black","grey")
             ,beside=F)
labs <- as.vector(colnames(counts))
text(cex=1.1, x=h, y=-17.5,labs, xpd=TRUE, srt=90)
h<-legend("topright", legend=rownames(counts), title="Alteration Type:", fill=c("black","grey"),cex=1.5)
dev.off()

