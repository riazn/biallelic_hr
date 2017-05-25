#clear window: 
rm(list=ls())
#Genes

gene_list_hr <- read.delim("Supplementary_Files/Master_List.txt", header=F)
gene_list_hr <- gsub(" ","",as.character(gene_list_hr$V1))
length(gene_list_hr)
Genes <- gene_list_hr

##OV,BRCA,PRAD with samples:
#LST >= 15
 

lst_brca <- read.delim("Supplementary_Files/LSTs/LST-BRCA.dat", header=F,sep="\t")
lst_ov <- read.delim("Supplementary_Files/LSTs/LST-OV.dat", header=F,sep="\t")
lst_prad <- read.delim("Supplementary_Files/LSTs/LST-PRAD.dat", header=F,sep="\t")

##High LST
lst_high <- c(as.character(lst_brca$V2[which(lst_brca$V3 >= 15)]),as.character(lst_ov$V2[which(lst_ov$V3 >=15)]),as.character(lst_prad$V2[which(lst_prad$V3>=15)]))
lst_high <- gsub(" ","",lst_high)

length(lst_high)

##OV,BRCA,PRAD with samples:
##Patients with dominant signature 3

###
Score <- read.delim("Supplementary_Files/Scores_Combined_BRCA_OV_PRAD_allLST_allBRCA.Signatures.txt",header=T,sep="\t")
dim(Score)

Mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(Mat)

samples <- Score[which(as.character(Score$Id)%in%gsub("\\.","-",colnames(Mat))),]
dim(samples)

maxim <- apply(samples[,3:15],1,function(x) x[which(x==max(x))])

samples$dominant <- NA
for(i in 1:nrow(samples)){
  
  if(length(maxim[[i]])==1){
    samples$dominant[i] <- names(maxim[[i]])}
}
dim(samples)
head(samples)

## 
samples_dom <- samples[which(samples$dominant%in%"Sig.3"),]
##Dominant MutSig3
m3_high <- gsub(" ","",samples_dom$Id)
length(m3_high)

## Select cases with HRD phenotype: LST >=15 OR Dominant MutSig3
## Take union ()

hrd <- union(lst_high,m3_high)
length(unique(hrd))


##Upload matrix

Mat <- read.delim("Supplementary_Files/Matrix_Biallelic_Monoallelic_Pathogenic_VUS_All_cancers_Mutation_Types_Paper.txt",header=T,sep="\t")
dim(Mat)


Mat_HBOCs_hrd <- Mat[,which(gsub("\\.","-",colnames(Mat))%in%hrd),]
dim(Mat_HBOCs_hrd)


## Cbio Oncoprinter
##Germline LOH HBOC cancers




#AMP --> Germ patho LOH mat==1
Dat=NULL
for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
 
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==1)])
  Alteration <- "AMP"
  Type <- "CNA"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  if(is.null(Dat)){Dat <- dat} else {Dat <- rbind(Dat,dat)}
}

dim(Dat)
#Germline Truncating: Add tag --> Trunc


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==1)])
  Alteration <- "p."
  Type <- "TRUNC"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  Dat <- rbind(Dat,dat)
  
}

dim(Dat)

##Somatic LOH mat == 5

for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==5)])
  Alteration <- "HOMDEL"
  Type <- "CNA"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  if(is.null(Dat)){Dat <- dat} else {Dat <- rbind(Dat,dat)}
}

dim(Dat)
#Somatic Truncating: Add tag --> Trunc


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==5)])
  Alteration <- "p."
  Type <- "TRUNC"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  Dat <- rbind(Dat,dat)
  
}

dim(Dat)
###Germline compund het mat==2


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==2)])
  Alteration <- "AMP"
  Type <- "CNA"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  if(is.null(Dat)){Dat <- dat} else {Dat <- rbind(Dat,dat)}
}

dim(Dat)
#Germline Truncating: Add tag --> Inframe


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==2)])
  Alteration <- "p."
  Type <- "INFRAME"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  Dat <- rbind(Dat,dat)
  
}

dim(Dat)

### Pathogenic somatic compund het: mat==6

for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==6)])
  Alteration <- "HOMDEL"
  Type <- "CNA"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  if(is.null(Dat)){Dat <- dat} else {Dat <- rbind(Dat,dat)}
}

dim(Dat)
#Somatic Truncating: Add tag --> Inframe


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==6)])
  Alteration <- "p."
  Type <- "INFRAME"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  Dat <- rbind(Dat,dat)
  
}

dim(Dat)


##VUS Somatic: somatic mat==7,8,3 

for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==7 | Mat_HBOCs_hrd[Gene,]==8 | Mat_HBOCs_hrd[Gene,]==3)])
  Alteration <- "HOMDEL"
  Type <- "CNA"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  if(is.null(Dat)){Dat <- dat} else {Dat <- rbind(Dat,dat)}
}

dim(Dat)
#Somatic VUS: Add tag --> Missense


for(k in 1:nrow(Mat_HBOCs_hrd))
{
  
  dat=NULL
  
  Gene <- as.character(rownames(Mat_HBOCs_hrd)[k])
  print(Gene)
  id <- gsub("\\.","-",colnames(Mat_HBOCs_hrd)[which(Mat_HBOCs_hrd[Gene,]==7 | Mat_HBOCs_hrd[Gene,]==8 | Mat_HBOCs_hrd[Gene,]==3)])
  Alteration <- "p."
  Type <- "MISSENSE"
  
  dat <- data.frame(Sample=id,Gene=rep(Gene,length(id)),Alteration=rep(Alteration,length(id)),Type=rep(Type,length(id)))
  Dat <- rbind(Dat,dat)
  
}

dim(Dat)
head(Dat)

##Number(bi allelic alt) >= 2 mat == 1,2,5,6

Mat <- Mat_HBOCs_hrd
dim(Mat)

for(k in 1:nrow(Mat)){
  Gene <- rownames(Mat)[k];
  print(Gene);
  Mat[Gene,which(Mat[Gene,]==2)] <- 1
  Mat[Gene,which(Mat[Gene,]==3)] <- 0
  Mat[Gene,which(Mat[Gene,]==4)] <- 0
  Mat[Gene,which(Mat[Gene,]==5)] <- 1
  Mat[Gene,which(Mat[Gene,]==6)] <- 1
  Mat[Gene,which(Mat[Gene,]==7)] <- 0
  Mat[Gene,which(Mat[Gene,]==8)] <- 0
  Mat[Gene,which(Mat[Gene,]==9)] <- 0
  Mat[Gene,which(Mat[Gene,]==10)] <- 0
}
dim(Mat)

length(which(Mat==1))
res <- apply(Mat,1,sum)
res <- data.frame(res)
cols <- rownames(res)[order(res$res,decreasing = T)]
res <- res[order(res$res,decreasing = T),]
res <- data.frame(gene=cols,count=res)
head(res)

dim(res[which(res$count >=2),])
top_16_bi_patho_genes <- res$gene[1:16]
top_16_bi_patho_genes

##Restric HBOC hrd matrix to those 16:
##Read Oncoprint and top 16 altered genes in HBOCs high LST>=15 OR Dominant MutSig3

Dat_top16 <- Dat[which(Dat$Gene%in%top_16_bi_patho_genes),]
dim(Dat_top16)
length(unique(Dat_top16$Gene))

write.table(Dat_top16,"Results_Figures_and_P_Values/Oncoprint_HBOCs_top16_LST15_OR_Dominant_MutSig3.txt", col.names=T,row.names=F,quote = F,sep="\t")

##Use cbio portal --> tools --> Oncoprinter:
##http://www.cbioportal.org/

