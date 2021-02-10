
##### #1: Summarise counts per gene ########## 

#source("https://bioconductor.org/biocLite.R")
#biocLite("tximportData")
#install.packages("readr")
library(tximportData)
library(readr)
library(rhdf5)

# read in pre-constructed tx2gene table (transcript to gene table)
tx2gene <- read.table("../philippa_resources/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt", header=T)
head(tx2gene)

files <- c('1.1','1.2','1.3','2.1','2.2','2.3','3.1','3.2','3.3')
names(files) <-c('1.1','1.2','1.3','2.1','2.2','2.3','3.1','3.2','3.3')

samples <- read.table(file.path(htseqDir,"samples.txt"), header=T)

library(tximport)
# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

samples <- read.table(file.path(htseqDir,"samples.txt"), header=T)


for(i in samples){
  # save counts summarised per gene
  write.table(txi$counts[,i], file=paste(i,"_count.tsv",sep=''),sep = "\t",quote=F,col.names = F)
}
# to see tpm summarised per gene
#head(txi$abundance)
#colnames(txi$abundance)

# save tpm summarised per gene
write.table(txi$abundance, file="_tpm.tsv",sep = "\t")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="_gene_lengths.csv")


