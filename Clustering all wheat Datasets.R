####################################################################
#Library
####################################################################
library(DESeq2)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(VennDiagram)
library(ggrepel)
library(rhdf5)
library(tximport)

####################################################################
# Set work environment and convert gene IDs to transcript IDs for txi import to DESeq2
####################################################################

setwd("~/Documents/Wheat/Transcriptomics_analysis")
htseqDir<-getwd()
samples <- read.table(file='All_datasets/Samples.txt', header=T)

##  read in pre-constructed tx2gene table (transcript to gene table)
##  If putting new cDNA into library for kallisto remember to do this in the transcripts to genes file too

tx2gene <- read.table("philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)
files <- as.character(samples$filename)
names(files) <- c("TE1","TE2","TE3","TM1","TM2","TM3","NE1","NE2","NE3","NIL7B_1","NIL7A_1","NIL7A_2","NIL7A_3",
                  "NIL1B_1","NIL1B_2","NIL1B_3","NIL1A_1","NIL7B_2","NIL7B_3","NIL1A_2","NIL1A_3",
                "Axillary_Meristem_1","Axillary_Meristem_2","Double_Ridge_Stage_1","Double_Ridge_Stage_2",
                  "Floral_Meristem_1","Floral_Meristem_2","Tetrads_Stage_1","Tetrads_Stage_2","Non_dividing_microspore_1",
                  "Non_dividing_microspore_2","Non_dividing_microspore_3","Dividing_microspore_1","Dividing_microspore_2","Dividing_microspore_3")

##  If changing pseudoaligner i.e. Salmon, change type=""
setwd("~/Documents/Wheat/Transcriptomics_analysis")
htseqDir<-getwd()
samples <- read.table(file='All_datasets/Samples_without_Svilena.txt', header=T)

##  read in pre-constructed tx2gene table (transcript to gene table)
##  If putting new cDNA into library for kallisto remember to do this in the transcripts to genes file too

tx2gene <- read.table("philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)
files <- as.character(samples$filename)
names(files) <- c("TE1","TE2","TE3","TM1","TM2","TM3","NE1","NE2","NE3","NIL7B_1","NIL7A_1","NIL7A_2","NIL7A_3",
                  "NIL1B_1","NIL1B_2","NIL1B_3","NIL1A_1","NIL7B_2","NIL7B_3","NIL1A_2","NIL1A_3",
                  "Axillary_Meristem_1","Axillary_Meristem_2","Double_Ridge_Stage_1","Double_Ridge_Stage_2",
                  "Floral_Meristem_1","Floral_Meristem_2","Tetrads_Stage_1","Tetrads_Stage_2")

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE,
                             countsFromAbundance=c("scaledTPM"))

head(txi.kallisto.tsv$abundances)


##  If changing pseudoaligner i.e. Salmon, change type=""
setwd("~/Documents/Wheat/Transcriptomics_analysis")
htseqDir<-getwd()
samples <- read.table(file='All_datasets/Samples_only_meristem_and_NIL.txt', header=T)

##  read in pre-constructed tx2gene table (transcript to gene table)
##  If putting new cDNA into library for kallisto remember to do this in the transcripts to genes file too

tx2gene <- read.table("philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)
files <- as.character(samples$filename)
names(files) <- c("NIL7B_1","NIL7A_1","NIL7A_2","NIL7A_3",
                  "NIL1B_1","NIL1B_2","NIL1B_3","NIL1A_1","NIL7B_2","NIL7B_3","NIL1A_2","NIL1A_3",
                  "Axillary_Meristem_1","Axillary_Meristem_2","Double_Ridge_Stage_1","Double_Ridge_Stage_2",
                  "Floral_Meristem_1","Floral_Meristem_2","Tetrads_Stage_1","Tetrads_Stage_2")

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE,
                             countsFromAbundance=c("scaledTPM"))

head(txi.kallisto.tsv$abundances)

####################################################################
# DESeq2 importing and normalising 
####################################################################

##  Read in the results from the analysis (the counts files)

ddsHTSeq<-DESeqDataSetFromTximport(txi=txi.kallisto.tsv,
                                   colData=samples,
                                   design =~condition)

##  design <- test everything in relations to condition
##  if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
##  Kallisto gives out TPM values but DESeq2 can only deal with integer count values so it will change them (or take the est. count column)

dds<-DESeq(ddsHTSeq)

filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.05,]
  return(res2)
}

##  alpha=p.value, lfcThreshold=log fold change, pAdjust..=hypothesis correction
##  only really appropriate to change lfcThreshold in my opinion, maybe pAdjust but p-value is just changing the rate of false positives

resultsNames(dds)
NIL7B_NIL7A_DEGs = results(dds, contrast= c("condition", "NIL7A", "NIL7B"), alpha = 0.05, pAdjustMethod = "BH",lfcThreshold = 0)
NIL7B_NIL7A_DEG = filter_degs(NIL7B_NIL7A_DEGs)

NIL1B_NIL1A_DEGs = results(dds, contrast= c("condition", "NIL1A", "NIL1B"), alpha = 0.05, pAdjustMethod = "BH",lfcThreshold = 0)
NIL1B_NIL1A_DEG = filter_degs(NIL1B_NIL1A_DEGs)

TM_NE_DEGs = results(dds, contrast= c("condition", "TM", "NE"), alpha = 0.05, pAdjustMethod = "BH",lfcThreshold = 0)
TM_NE_DEG = filter_degs(TM_NE_DEGs)

Dividing_DEGs = results(dds, contrast= c("condition", "Dividing_microspore", "Non_dividing_microspore"), alpha = 0.05, pAdjustMethod = "BH",lfcThreshold = 0)
Dividing_DEG = filter_degs(Dividing_DEGs)
###################################################################
#Preparing lists of diffrentially expressed genes for online tools
###################################################################

################## All DEGs ##################

# Extracting Gene IDs

NIL7B_NIL7A_DEG_ID <- rownames(NIL7B_NIL7A_DEG)
NIL1B_NIL1A_DEG_ID <- rownames(NIL1B_NIL1A_DEG)

################## Upregulated DEGs ##################

# Extracting upregulated genes

NIL7B_NIL7A_upregulated_DEG<-NIL7B_NIL7A_DEG[NIL7B_NIL7A_DEG$log2FoldChange>0,]
NIL1B_NIL1A_upregulated_DEG<-NIL1B_NIL1A_DEG[NIL1B_NIL1A_DEG$log2FoldChange>0,]
TM_NE_upregulated_DEG<-TM_NE_DEG[TM_NE_DEG$log2FoldChange>0,]
Dividing_upregulated_DEG<-Dividing_DEG[Dividing_DEG$log2FoldChange>0,]

##  Gene IDs

NIL7B_NIL7A_upregulated_DEG_ID <- rownames(NIL7B_NIL7A_upregulated_DEG)
NIL1B_NIL1A_upregulated_DEG_ID <- rownames(NIL1B_NIL1A_upregulated_DEG)

################## Downregulated DEGs ##################

# Extracting downregulated genes

NIL7B_NIL7A_downregulated_DEG<-NIL7B_NIL7A_DEG[NIL7B_NIL7A_DEG$log2FoldChange<0,]
NIL1B_NIL1A_downregulated_DEG<-NIL1B_NIL1A_DEG[NIL1B_NIL1A_DEG$log2FoldChange<0,]


##  Gene IDs

NIL7B_NIL7A_downregulated_DEG_ID <- rownames(NIL7B_NIL7A_downregulated_DEG)
NIL1B_NIL1A_downregulated_DEG_ID <- rownames(NIL1B_NIL1A_downregulated_DEG)

################## Shared DEGs Heatmap ##################

## Find shared DEGs

Shared_DEGs_IDs<- NIL7B_NIL7A_DEG_ID[NIL7B_NIL7A_DEG_ID %in% NIL1B_NIL1A_DEG_ID]
Shared_upregulated_DEGs_IDs<- NIL7B_NIL7A_upregulated_DEG_ID[NIL7B_NIL7A_upregulated_DEG_ID %in% NIL7B_NIL7A_upregulated_DEG_ID]
Shared_downregulated_DEGs_IDs<- NIL1B_NIL1A_downregulated_DEG_ID[NIL1B_NIL1A_downregulated_DEG_ID %in% NIL1B_NIL1A_downregulated_DEG_ID]

library(pheatmap)

counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
Shared_counts <- counts[Shared_DEGs_IDs,]
heatmap <- pheatmap((log2(Shared_counts+1)), scale = "row",border_color=NA,show_rownames = F,
         show_colnames=T,main = 'Shared DEGs expression across samples')


plot_tree(heatmap)
produce_clusters(heatmap,7.8)
cluster<-1
cluster_to_heatmap(cluster)

specific_genes <- unique((melt_counts[melt_counts$cluster==cluster,])[,2])
specific_genes

write.table(specific_genes,'cluster_1_specific_genes.txt',row.names = F,quote=F,col.names = F)


counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
Shared_counts <- counts[Shared_DEGs_IDs,]
heatmap <- pheatmap((log2(Shared_counts[,1:29]+1)), scale = "row",border_color=NA,show_rownames = F,
                    show_colnames=T,main = 'Shared DEGs expression across samples')


plot_tree(heatmap)
produce_clusters(heatmap,10.2)
cluster<-15
cluster_to_heatmap(cluster)

TM_vs_NE_DEG_IDs <- read.delim("~/Documents/Wheat/Transcriptomics_analysis/RKD/analysis/DEGS/Upregulated/TM_vs_NE_upregulated_DEGs_IDs_txt")
Specified <- as.data.frame(specific_genes)
TM_Spike_NIL_Cluster1 <- TM_vs_NE_DEG_IDs[TM_vs_NE_DEG_IDs %in% Specified]

write.table(specific_genes,'cluster_15_specific_genes.txt',row.names = F,quote=F,col.names = F)

vp = venn.diagram(list(TM_vs_NE_DEG_IDs = row.names(TM_vs_NE_DEG_IDs),
                       Specified = row.names(Specified)),
                  fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cex = 2, lty =2, 
                  filename = NULL)

grid.newpage()
grid.draw(vp)
###################################################################
#TF terms
###################################################################
TFs<- read.delim("~/Documents/R scripts and projects/Projects/RKD_grasses/philippa_resources/TFs/All.txt")

###################################################################
#Gage and list construction
###################################################################
library(gage)

#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(NIL1B_NIL1A_DEGs$log2FoldChange))

###################################################################
list <- list()
for(i in 1:58){
  
  TF_class <- as.character(unique(TFs$superfamily))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name),TFs$superfamily),2]
}
names(list)<-TF_class[1:58]

#Run GAGE command for all leaky and induced expressed transgenics

Enriched_NILA_vs_NILB <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(10,14,15,16,18,19),samp=c(11,12,13,17,20,21),
                 rank.test = T,me.dir = F,
                 set.size=c(1,800), compare="paired",same.dir = T)

Enriched_Spike_vs_NILB <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(10,14,15,16,18,19),samp=c(22:29),
                              rank.test = T,me.dir = F,
                              set.size=c(1,800), compare="unpaired",same.dir = T)

Enriched_Dividing_Microspore <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(30:32),samp=c(33:35),
                               rank.test = T,me.dir = F,
                               set.size=c(1,800), compare="paired",same.dir = T)
###################################################################
#corrplot
###################################################################
RKD_wheat_TF_enrichment <- read.delim('~/Documents/Wheat/Transcriptomics_analysis/RKD/analysis/RKD_TF_gene_enrichment.txt')

RKD_wheat <- RKD_wheat_TF_enrichment[ order(row.names(RKD_wheat_TF_enrichment)), ]
RKD_wheat<-RKD_wheat[,c(1:4)]

Enriched_NILA_vs_NILB <- Enriched_NILA_vs_NILB$stats[ order(row.names(Enriched_NILA_vs_NILB$stats)), ]

Enriched_Spike_vs_NILB <- Enriched_Spike_vs_NILB$stats[ order(row.names(Enriched_Spike_vs_NILB$stats)), ]

Enriched_Dividing_Microspore <- Enriched_Dividing_Microspore$stats[ order(row.names(Enriched_Dividing_Microspore$stats)), ]

Tillering_Spike <- cbind(Enriched_NILA_vs_NILB, Enriched_Spike_vs_NILB, RKD_wheat,Enriched_Dividing_Microspore)[,c(1,8,17,21)]
Tillering_Spike <-Tillering_Spike[-36,]

colnames(Tillering_Spike) <- c("NIL","Spike development", "RKD","Microspore")

pheatmap(Tillering_Spike)

###########################################
#GO term analysis 
############################################

TM_DEGs_GO_term <- read.csv('All_datasets/Cluster_1_inflorescence_NIL.csv')
TM_DEGs_GO_term <- TM_DEGs_GO_term[,c(1:4)]
TM_DEGs_GO_term$Ratio <- TM_DEGs_GO_term$Genes.in.list/TM_DEGs_GO_term$Total.genes
as.factor(TM_DEGs_GO_term$Enrichment.FDR)
TM_DEGs_GO_term$Enrichment.FDR <- -log10(TM_DEGs_GO_term$Enrichment.FDR)


ggplot(TM_DEGs_GO_term, mapping=aes(Functional.Category,Ratio,colour=Enrichment.FDR)) +
  coord_flip()+
  geom_point(aes(size = Genes.in.list)) +
  scale_colour_gradient(high="blue",low="red")




