#clear window: 
rm(list=ls())


## Per gene accross cancers mut types

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)

## Load matrix
mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat)


List_All_Genes_HR=NULL

  
  for(k in 1:length(gene_list_hr)){
    
    
    gene <- gene_list_hr[k];
    
    print(gene);
    
    germ_bi_patho <- 0
    germ_bi_patho_double <- 0
    germ_bi_vus <- 0
    
    som_bi_patho <- 0
    som_bi_patho_double <- 0
    som_bi_vus <- 0
    
    som_loh <- 0
    
    
    germ_bi_patho_num <- length(which(mat[gene,]==1));
    germ_bi_patho_double_num <- length(which(mat[gene,]==2));
    germ_bi_vus_num <- length(which(mat[gene,]==3));
    
    som_bi_patho_num <- length(which(mat[gene,]==5));
    som_bi_patho_double_num <- length(which(mat[gene,]==6));
    som_bi_vus_num <- length(which(mat[gene,]==7));
    
    som_loh_num <- length(which(mat[gene,]==8));
    
    germ_bi_patho <- germ_bi_patho + germ_bi_patho_num;
    germ_bi_patho_double <- germ_bi_patho_double + germ_bi_patho_double_num;
    germ_bi_vus <- germ_bi_vus + germ_bi_vus_num;
    
    som_bi_patho <- som_bi_patho + som_bi_patho_num;
    som_bi_patho_double <- som_bi_patho_double + som_bi_patho_double_num;
    som_bi_vus <- som_bi_vus + som_bi_vus_num;
    
    som_loh <- som_loh + som_loh_num
  
  res <- data.frame(germ_biallelic_pathogenic=germ_bi_patho,germ_biallelic_pathogenic_double=germ_bi_patho_double,germ_biallelic_vus=germ_bi_vus,som_biallelic_pathogenic=som_bi_patho,som_biallelic_pathogenic_double=som_bi_patho_double,som_biallelic_vus=som_bi_vus,som_loh=som_loh)
  if(is.null(List_All_Genes_HR)) {List_All_Genes_HR <- res} else {
    List_All_Genes_HR <- rbind(List_All_Genes_HR,res)}
}

rownames(List_All_Genes_HR) <- gene_list_hr

res <- List_All_Genes_HR
colnames(res) <- c("Germline Trunc LOH","Germline Trunc Somatic Trunc","Germline Trunc Somatic not Trunc","Somatic Trunc LOH","Somatic Trunc Somatic Trunc","Somatic Trun Somatic not Trunc","Somatic not Trunc LOH")
head(res)
dim(res)

##Count pathogenic events
Count=NULL
for(k in 1:dim(res)[1])
{
  gene <- rownames(res)[k]
  print(gene)
  count <- res[gene,"Germline Trunc LOH"]+res[gene,"Germline Trunc Somatic Trunc"]+res[gene,"Somatic Trunc LOH"]+res[gene,"Somatic Trunc Somatic Trunc"]
  Count <- c(Count,count)
}

count <- data.frame(Count)
rownames(count) <- rownames(res)
res_count <- cbind(res,count)
res_count <- data.frame(res_count)
res_count <- res_count[order(res_count$Count,decreasing = T),]
head(res_count)
res_count[1:25,]
##Merge compund hets, germ and som

res_merged <- data.frame(res_count$Germline.Trunc.LOH,res_count$Germline.Trunc.Somatic.Trunc + res_count$Germline.Trunc.Somatic.not.Trunc,res_count$Somatic.Trunc.LOH,res_count$Somatic.Trunc.Somatic.Trunc + res_count$Somatic.Trun.Somatic.not.Trunc,res_count$Somatic.not.Trunc.LOH,res_count$Count)
head(res_merged)
rownames(res_merged) <- rownames(res_count)


table <- data.matrix(t(res_merged))
head(table)
rownames(table) <-  c("Germ Patho LOH","Germ Comp Het","Som Patho LOH","Som Comp Het","Som VUS LOH","Count")
head(table)

#Select Top 25

table <- table[,1:25]
rownames(table) <-  c("Germ Patho LOH","Germ Comp Het","Som Patho LOH","Som Comp Het","Som VUS LOH","Count")

#Plot

counts <- table[-6,]
pdf('Results_Figures_and_P_Values/Supplementary_Fig1c.pdf')
h <- barplot(counts,xaxt="n", main="",cex.main=0.9, ylab="Altered Tumors", cex.lab=1.5
             ,col=c("black","blue","grey","lightblue","darkgreen")
             ,beside=F)
labs <- as.vector(colnames(counts))
text(cex=1.3, x=h, y=-11,labs, xpd=TRUE, srt=90)
h<-legend("topright", legend=rownames(counts), title="Alteration Type:", fill=c("black","blue","grey","lightblue","darkgreen"),cex=1)
dev.off()

