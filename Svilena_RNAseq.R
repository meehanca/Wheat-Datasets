library(DESeq2)
htseqDir<-getwd()

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("./philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)

files <- c('./raw_svilena_counts/1.1','./raw_svilena_counts/1.2','./raw_svilena_counts/1.3',
           './raw_svilena_counts/2.1','./raw_svilena_counts/2.2','./raw_svilena_counts/2.3',
           './raw_svilena_counts/3.1','./raw_svilena_counts/3.2','./raw_svilena_counts/3.3')
names(files) <-c('1.1','1.2','1.3','2.1','2.2','2.3','3.1','3.2','3.3')

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

samples <- read.table(file='samples.txt', header=T)
##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromTximport(txi=txi,
                                   colData=samples,
                                   design =~condition)
## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

plotDispEsts(dds, main="Dispersion plot")
####################################################################
# Do PCA
####################################################################
#principal component analysis

(rowVars(assay(rld)))
vst = vst(ddsHTSeq)

v <- plotPCA(vst, intgroup=c("condition"))
v<- v+ geom_label_repel(aes(label = name))
v
pcaData <- DESeq2::plotPCA(vst, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA_parents.pdf", height = 6, width = 6)
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  #scale_colour_manual(name="",values = c("a12"="goldenrod2", "gd33"="darkslateblue", "f1"="saddlebrown"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme_bw()
#dev.off()

####################################################################
#Plotting Reps
####################################################################
plot_reps =  function(dds,x=1,y=2,cond_choice=1, cond='condition'){
  ##  Estimate the size factors for normalisation
  dds<-estimateSizeFactors(dds)
  
  ## Extract the normalised counts for the condition you want
  rep_values<- counts(dds, normalized=TRUE)[,dds[[cond]]==cond_choice]
  
  # Take logs of these values
  vals <- log2(rep_values[,c(x,y)] + 0.5)
  # And plot
  plot(vals,pch=16, cex=0.4,xlab=paste('rep',x),ylab=paste('rep',y))
  grid(col = "darkgray", lty = "solid",lwd = par("lwd"), equilogs = TRUE)
  title(paste("Comparison of",cond_choice,"replicates"))
}

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="untreated")
plot_reps(dds, x=1, y=3, cond_choice="untreated")
plot_reps(dds, x=2, y=3, cond_choice="untreated")

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="pretreatment")
plot_reps(dds, x=1, y=3, cond_choice="pretreatment")
plot_reps(dds, x=2, y=3, cond_choice="pretreatment")

par(mfrow = c(3,1))

plot_reps(dds, x=1, y=2, cond_choice="treatment")
plot_reps(dds, x=1, y=3, cond_choice="treatment")
plot_reps(dds, x=2, y=3, cond_choice="treatment")

####################################################################
#DEGs
####################################################################
filter_degs <- function(res){
  summary(res)
  res2 = res[!(is.na(res$padj)),]
  res2 = res2[res2$padj < 0.01,]
  return(res2)
}

untreated_vs_treated_DEGs = results(dds, contrast= c("condition", "untreated", "treatment"), alpha = 0.01, pAdjustMethod = "BH",lfcThreshold=1)
untreated_vs_treated_DEG = filter_degs(untreated_vs_treated_DEGs)

pretreated_vs_treated_DEGs = results(dds, contrast= c("condition", "pretreatment", "treatment"), alpha = 0.01, pAdjustMethod = "BH",lfcThreshold=1)
pretreated_vs_treated_DEG = filter_degs(pretreated_vs_treated_DEGs)

summary(untreated_vs_pretreated_DEG)
head(untreated_vs_pretreated_DEG)

summary(untreated_vs_treated_DEG)
head(untreated_vs_treated_DEG)

summary(pretreated_vs_treated_DEG)
head(pretreated_vs_treated_DEG)
####################################################################
#MA Plots
####################################################################
par(mfrow=c(3,1))
DESeq2::plotMA(untreated_vs_pretreated_DEGs, ylim=c(-10,15), main='untreated_vs_pretreated_DEGs')
DESeq2::plotMA(untreated_vs_treated_DEGs, ylim=c(-10,15), main='untreated_vs_treated_DEGs')
DESeq2::plotMA(pretreated_vs_treated_DEGs, ylim=c(-10,15), main='pretreated_vs_treated_DEGs')
####################################################################
#Volcano
####################################################################
log10.pval <- -log10(untreated_vs_pre-treated_DEGs$padj)
log2.fc    <- untreated_vs_pre-treated_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,75),
     main='untreated_vs_pre-treated_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)

log10.pval <- -log10(untreated_vs_treated_DEGs$padj)
log2.fc    <- untreated_vs_treated_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,75),
     main='untreated_vs_treated_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)

log10.pval <- -log10(pre-treated_vs_treated_DEGs$padj)
log2.fc    <- pre-treated_vs_treated_DEGs$log2FoldChange
plot(log2.fc,log10.pval,
     xlab="log2 (fold change)",
     ylab="-log10 (p-value)",
     col="black",
     xlim=c(-10,10),
     ylim=c(0,75),
     main='pre-treated_vs_treated_DEGs')
abline(h = -log10(0.05),col='red',lwd=1.5)
abline(v=-log2(2),col='blue',lwd=1.5)
abline(v=log2(2),col='blue',lwd=1.5)
####################################################################
#Venn
####################################################################
vp = venn.diagram(list(untreated_vs_pretreated = row.names(untreated_vs_pretreated_DEG),
                       untreated_vs_treated = row.names(untreated_vs_treated_DEG)),
                  fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cex = 2, lty =2, 
                  filename = NULL)

grid.newpage()
grid.draw(vp)

vp = venn.diagram(list(pretreated_vs_treated = row.names(pretreated_vs_treated_DEG),
                       untreated_vs_treated = row.names(untreated_vs_treated_DEG)),
                  fill = c("red", "blue"),
                  alpha = c(0.5, 0.5), cex = 2, lty =2, 
                  filename = NULL)

grid.newpage()
grid.draw(vp)

####################################################################
#Heatmap
####################################################################
###Only provided with lines to make a heatmap for VIM1 mutant, try to figure out how to make heatmaps for the other 2
counts = counts(dds , normalized = TRUE)
counts <-  counts[apply(counts, MARGIN = 1, FUN = function(x) sd(x) != 0 ),]#it removes genes that are not express and have no variance
colnames(counts) <- c("untreated1","untreated2","untreated3","pretreated1","pretreated2","pretreated3", "treated1", "treated2", "treated")
pretreated_vs_treated_DEG_counts <- (counts[rownames(pretreated_vs_treated_DEG),])
pheatmap((log2(pretreated_vs_treated_DEG_counts+1)), scale = "row",border_color=NA,show_rownames = F, main = 'pretreated vs treated DEGs expression across samples')

#we are taking log2 to reduce the difference, srink the values 
# you can remove the scale if you want to remove that scale you will see the massive difference 
# +1 is important as well, but I don't know why :-(
#cluster_rows = F
#cluster_cols = F
# these arguments are by default set as T, if you change them to false interesting things can happen 

### try a random forest walk? 
#did you bring your boots?

require(randomForest)

# there are some random features to randomForest
# we can use a 'random seed' to make the random features reproducible
# you may use any seed you like, but for this example let's use the same one
set.seed(1234)

# lets try making a model
# we will predict FPKM (x)
# using genotype (y)
# and a small number of decision trees (ntree)

samples <- read.table(file='sampleTable2.txt', header=T)
##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromTximport(txi=txi,
                                   colData=samples,
                                   design =~condition)
## design<- you say to the test to do everything in relations to condition
## if you have more than one conditions you want to differentiate (for example different genotypes) you change design = ~  condition + genotype
##  And perform the analysis (details in the manual)

dds<-DESeq(ddsHTSeq)

# x
t(counts(dds))

# y
colData(dds)$condition

# let's calculate the relative importance of all genes involved in making the decision tree

rf.50 <- randomForest(x=t(counts(dds)), y=colData(dds)$condition, ntree=11000, importance=T)
print(rf.50)

# plot importances

varImpPlot(rf.50, n.var=50, main="Top 25 Influential Transcripts")

# now to refine the number of genes considered in each tree
# the randomForest package comes with a built in function that will try to optimize
# by trying mtry values higher or lower than the starting one, and looking for improvement
# when it can't improve by calculating a new mtry tree, it stops
mtry <- tuneRF(x = t(counts(dds)), y = colData(dds)$condition, ntreeTry=50, mtryStart=100,
               stepFactor=5,improve=1, trace=TRUE, plot=TRUE) #increases/decreases by step factor and checks for improve %improvement
print(mtry)

#extract influence metrics from random forest generation
influence <- importance(rf.50)
#order by error and output
head(influence[order(influence[,5], decreasing=T),])
influence.order <- influence[order(influence[,5], decreasing=T),]
write.table(influence.order, 
            file="randomForest_importance.13S.5C.txt", sep="\t", quote=F)

plotCounts(dds, gene="TraesCS3D02G310600", intgroup="condition")
plotCounts(dds, gene="TraesCS2A02G218900", intgroup="condition") 
plotCounts(dds, gene="TraesCS2A02G414300", intgroup="condition") 
plotCounts(dds, gene="TraesCS5A02G314000", intgroup="condition")

###################################################################
#transcription factors
###################################################################

####Upregulated

TFs <- read.table('./TFs/shared_upregulated_TFs.txt',sep='\t',header=F)[,c(2,1)]
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 33
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`))) + ggtitle("Shared upregulated DEGs") +
  geom_bar(stat="count", position = "dodge",size=2) + theme_minimal() + theme(
    axis.text.x=element_text(angle = -90, hjust = 0),
    axis.ticks.x=element_blank()) + 
  scale_fill_manual(values =c(rep("grey",15))) + xlab(element_blank()) + ylab("Number of Transcription Factors")

####All TFs

TFs <- read.table('TF_shared.txt',sep=',',header=F)
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 33
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`),fill = TFs$`TF class`)) +
  geom_bar(stat="count", position = "dodge",size=2) + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
  scale_fill_manual(values =mycolors) + xlab("TF class") + ylab("Number of shared TFs of that class")

####Downregulated


TFs <- read.table('TFs_downregulated_shared_DEG.txt',sep=',',header=F)
colnames(TFs) <- c("1","TF class") 
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 33
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

ggplot(TFs, aes(factor(TFs$`TF class`),fill = TFs$`TF class`)) +
  geom_bar(stat="count", position = "dodge",size=2) + theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
  scale_fill_manual(values =mycolors) + xlab("TF class") + ylab("Number of downregulated TFs of that class")

###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

untreated_vs_pretreated_counts <- rownames(untreated_vs_pretreated_DEG)
untreated_vs_treated_counts <- rownames(untreated_vs_treated_DEG)
pretreated_vs_treated_counts <- rownames(pretreated_vs_treated_DEG)

write.table(untreated_vs_pretreated_counts, file = "untreated_vs_pretreated_counts.txt", quote = F, row.names = F, col.names = F)
write.table(untreated_vs_treated_counts, file = "untreated_vs_treated_counts.txt", quote = F, row.names = F, col.names = F)
write.table(pretreated_vs_treated_counts, file = "pretreated_vs_treated_counts.txt", quote = F, row.names = F, col.names = F)

#counts
untreated_vs_pretreated_DEG_counts <- cbind(rownames(untreated_vs_pretreated_DEG),as.data.frame(untreated_vs_pretreated_DEG))
untreated_vs_treated_DEG_counts <- cbind(rownames(untreated_vs_treated_DEG),as.data.frame(untreated_vs_treated_DEG))
pretreated_vs_treated_DEG_counts <- cbind(rownames(pretreated_vs_treated_DEG),as.data.frame(pretreated_vs_treated_DEG))

write.table(untreated_vs_pretreated_DEG_counts, file = "untreated_vs_pretreated_DEG_counts.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_treated_DEG_counts, file = "untreated_vs_treated_DEG_counts.txt", quote = F, row.names = T, col.names = T)
write.table(pretreated_vs_treated_DEG_counts, file = "pretreated_vs_treated_DEG_counts.txt", quote = F, row.names = T, col.names = T)

pretreated_vs_treated_upregulated_DEGs<-pretreated_vs_treated_DEG[pretreated_vs_treated_DEG$log2FoldChange>0,]
pretreated_vs_treated_downregulated_DEGs<-pretreated_vs_treated_DEG[pretreated_vs_treated_DEG$log2FoldChange<0,]

untreated_vs_treated_upregulated_DEGs<-untreated_vs_treated_DEG[untreated_vs_treated_DEG$log2FoldChange>0,]
untreated_vs_treated_downregulated_DEGs<-untreated_vs_treated_DEG[untreated_vs_treated_DEG$log2FoldChange<0,]

write.table(pretreated_vs_treated_upregulated_DEGs, file = "pretreated_vs_treated_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(pretreated_vs_treated_downregulated_DEGs, file = "pretreated_vs_treated_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_treated_upregulated_DEGs, file = "untreated_vs_treated_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_treated_downregulated_DEGs, file = "untreated_vs_treated_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

Shared_DEGs<-subset(untreated_vs_treated_DEG, rownames(untreated_vs_treated_DEG) %in% rownames(pretreated_vs_treated_DEG))

Shared_up_DEGs<-Shared_DEGs[Shared_DEGs$log2FoldChange>0,]
Shared_down_DEGs<-Shared_DEGs[Shared_DEGs$log2FoldChange<0,]

###################################################################
#TF terms
###################################################################
TFs<- read.delim("~/Documents/R scripts and projects/Projects/RKD_grasses/philippa_resources/TFs/All.txt")

###################################################################
#Gage and list construction
###################################################################
library(gage)
library(ggplot2)
#Exclude lowly expressed genes for GSEA
DESeq2_negative_gene_IDs <- is.na(as.data.frame(untreated_vs_treated_DEGs$log2FoldChange))

###################################################################
list <- list()
for(i in 1:58){
  
  TF_class <- as.character(unique(TFs$superfamily))
  TF_class_name <- TF_class[i]
  
  list[[i]] <- TFs[grep(paste(TF_class_name),TFs$superfamily),2]
}
names(list)<-TF_class[1:58]

#Run GAGE command for all leaky and induced expressed transgenics

Enriched <- gage(counts(dds)[!DESeq2_negative_gene_IDs,],list,ref=c(1:6),samp=c(7:9),
                 rank.test = T,me.dir = F,
                 set.size=c(1,800), compare="unpaired")

Enriched_greater <- Enriched$greater[1:18,1:5]
Enriched_lesser <- Enriched$less[1:14,1:5]

Enriched_write <- rbind(Enriched_greater,Enriched_lesser)
write.table(Enriched_write,file = '../../analysis/TF_gene_enrichment.txt', quote =F, sep= "\t")
write.table(Enriched_greater,file = '../../analysis/TF_greater_gene_enrichment.txt', quote =F, sep= "\t")

q.val <- -log10(Enriched_write[,4])
q.val[16:32] <- q.val[16:32]*-1

data<-data.frame(rownames(Enriched_write), q.val)
colnames(data) <- c("TF","q.val")

colours <- c(rep("indianred1",15), rep("royalblue",17))

library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

library(wesanderson)

ggplot(data, aes(TF, q.val, fill=q.val)) + geom_bar(stat="identity") +
  scale_fill_continuous(low="blue", high="red") +
  coord_flip()

ggplot(userData, aes(month, count, fill = count)) +
  geom_bar(stat = "identity") +
  scale_x_date() + 
  scale_fill_continuous(low="blue", high="red") +
  labs(x= "Time", y="Count")




# normalize and rlog-transform the data to get rid of most of the technical/uninteresting variation
rlog_ds <- rlog(ddsHTSeq)
# calculate variance (or whatever measure you prefer) per gene
rv <- rowVars(assay(rlog_ds))
# sort, so that the most variable genes will be on top of the object
rv <- order(rv, decreasing = TRUE)

