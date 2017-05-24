# Clear window: 
rm(list=ls())

##Generate data

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)

## load matrix

mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(mat)


List_All_Cancers_HR=NULL

  
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
  if(is.null(List_All_Cancers_HR)) {List_All_Cancers_HR <- res} else {
    List_All_Cancers_HR <- rbind(List_All_Cancers_HR,res)}
}

rownames(List_All_Cancers_HR) <- gene_list_hr
head(List_All_Cancers_HR)
write.table(List_All_Cancers_HR,"Supplementary_Files/Alterations_per_Gene_across_cancers_HR_Pathway_ordered_Pathogenic_Master_Paper.txt", col.names=T,row.names=T,quote=F,sep="\t")


#Figure 1A) Pathogenic Alterations per gene accross cancers: Germline black Somatic grey top 25 genes ordered by pathogenicity

patho_germ_som <- read.csv("Supplementary_Files/Alterations_per_Gene_across_cancers_HR_Pathway_ordered_Pathogenic_Master_Paper.txt", header=T,sep="\t")
res <- patho_germ_som
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

table <- data.matrix(t(res_count))


#Select Top 25

table <- table[,1:25]
rownames(table) <-  c("Germline Pathogenic with LOH","Germline Pathogenic with Somatic Pathogenic","Germline Pathogenic with Somatic VUS","Somatic Pathogenic with LOH","Somatic Pathogenic with Somatic Pathogenic","Somatic Pathogenic with Somatic VUS","Somatic VUS with LOH","Count")
##Merge germline and somatic

table_germ <- apply(table[1:2,],2,sum)
table_som <-  apply(table[4:5,],2,sum)
table_germ_som <- rbind(table_germ,table_som)
rownames(table_germ_som) <- c("Germline Biallelic Pathogenic","Somatic Biallelic Pathogenic")

#Plot
counts <- table_germ_som
pdf('Results_Figures_and_P_Values/Fig1a.pdf')
h <- barplot(counts,xaxt="n", main="",cex.main=0.9, ylab="Altered Tumors", cex.lab=1.5
             ,col=c("black","grey")
             ,beside=F)
labs <- as.vector(colnames(counts))
text(cex=1.2, x=h, y=-9.3,labs, xpd=TRUE, srt=90)
h<-legend("topright", legend=rownames(counts), title="Alteration Type:", fill=c("black","grey"),cex=1.5)
dev.off()

