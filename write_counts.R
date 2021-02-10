###################################################################
#Preparing lists of diffrentially expressed genes for online tools
####################################################################

################## All DEGs ##################

# Extracting Gene IDs

untreated_vs_treated_DEG_ID <- rownames(untreated_vs_treated_DEG)
untreated_vs_pretreated_DEG_ID <- rownames(untreated_vs_pretreated_DEG)
pretreated_vs_treated_DEG_ID <- rownames(pretreated_vs_treated_DEG)

# Writing tables - All DEG Counts

write.table(untreated_vs_treated_DEG, file = "DEGs/DEGs/untreated_vs_treated_DEG_counts.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_pretreated_DEG, file = "DEGs/DEGs/untreated_vs_pretreated_DEG_counts.txt", quote = F, row.names = T, col.names = T)
write.table(pretreated_vs_treated_DEG, file = "DEGs/DEGs/pretreated_vs_treated_DEG_counts.txt", quote = F, row.names = T, col.names = T)

# Writing tables - All DEG gene IDs

write.table(untreated_vs_treated_DEG_ID, file = "DEGS/DEGs/untreated_vs_treated_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(untreated_vs_pretreated_DEG_ID, file = "DEGs/DEGs/untreated_vs_pretreated_DEG_ID.txt", quote = F, row.names = F, col.names = F)
write.table(pretreated_vs_treated_DEG_ID, file = "DEGs/DEGs/pretreated_vs_treated_DEG_ID.txt", quote = F, row.names = F, col.names = F)

################## Upregulated DEGs ##################

# Extracting upregulated genes

untreated_vs_treated_upregulated_DEG<-untreated_vs_treated_DEG[untreated_vs_treated_DEG$log2FoldChange>0,]
pretreated_vs_treated_upregulated_DEG<-pretreated_vs_treated_DEG[pretreated_vs_treated_DEG$log2FoldChange>0,]
untreated_vs_pretreated_upregulated_DEG<-untreated_vs_pretreated_DEG[untreated_vs_pretreated_DEG$log2FoldChange>0,]

# Writing tables - Upregulated DEG counts

write.table(untreated_vs_treated_upregulated_DEG, file = "DEGs/Upregulated/untreated_vs_treated_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(pretreated_vs_treated_upregulated_DEG, file = "DEGs/Upregulated/pretreated_vs_treated_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_pretreated_upregulated_DEG, file = "DEGs/Upregulated/untreated_vs_pretreated_upregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

untreated_vs_treated_upregulated_DEG_ID <- rownames(untreated_vs_treated_upregulated_DEG)
untreated_vs_pretreated_upregulated_DEG_ID <- rownames(untreated_vs_pretreated_upregulated_DEG)
pretreated_vs_treated_upregulated_DEG_ID <- rownames(pretreated_vs_treated_upregulated_DEG)

# Writing tables - Upregulated DEG IDs

write.table(untreated_vs_treated_upregulated_DEG_ID, file = "DEGs/Upregulated/untreated_vs_treated_upregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(pretreated_vs_treated_upregulated_DEG_ID, file = "DEGs/Upregulated/pretreated_vs_treated_upregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(untreated_vs_pretreated_upregulated_DEG_ID, file = "DEGs/Upregulated/untreated_vs_pretreated_upregulated_DEGs_IDs_txt", quote = F, row.names = F, col.names = F)

################## Downregulated DEGs ##################

# Extracting downregulated genes

untreated_vs_treated_downregulated_DEG<-untreated_vs_treated_DEG[untreated_vs_treated_DEG$log2FoldChange<0,]
pretreated_vs_treated_downregulated_DEG<-pretreated_vs_treated_DEG[pretreated_vs_treated_DEG$log2FoldChange<0,]
untreated_vs_pretreated_downregulated_DEG<-untreated_vs_pretreated_DEG[untreated_vs_pretreated_DEG$log2FoldChange<0,]

# Writing tables - Downregulated DEG counts

write.table(untreated_vs_treated_downregulated_DEG, file = "DEGs/Downregulated/untreated_vs_treated_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(pretreated_vs_treated_downregulated_DEG, file = "DEGs/Downregulated/pretreated_vs_treated_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)
write.table(untreated_vs_pretreated_downregulated_DEG, file = "DEGs/Downregulated/untreated_vs_pretreated_downregulated_DEGs.txt", quote = F, row.names = T, col.names = T)

##  Gene IDs

untreated_vs_treated_downregulated_DEG_ID <- rownames(untreated_vs_treated_downregulated_DEG)
untreated_vs_pretreated_downregulated_DEG_ID <- rownames(untreated_vs_pretreated_downregulated_DEG)
pretreated_vs_treated_downregulated_DEG_ID <- rownames(pretreated_vs_treated_downregulated_DEG)

# Writing tables - Downregulated DEG counts

write.table(untreated_vs_treated_downregulated_DEG_ID, file = "DEGs/Downregulated/untreated_vs_treated_downregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(pretreated_vs_treated_downregulated_DEG_ID, file = "DEGs/Downregulated/pretreated_vs_treated_downregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(untreated_vs_pretreated_downregulated_DEG_ID, file = "DEGs/Downregulated/untreated_vs_pretreated_downregulated_DEGs_IDs_txt", quote = F, row.names = F, col.names = F)

################## Shared DEGs between (TE vs NE) and (TM vs NE) ##################

## Find shared DEGs

Shared_DEGs_IDs<- pretreated_vs_treated_DEG_ID[pretreated_vs_treated_DEG_ID %in% untreated_vs_treated_DEG_ID]
Shared_upregulated_DEGs_IDs<- pretreated_vs_treated_upregulated_DEG_ID[pretreated_vs_treated_upregulated_DEG_ID %in% untreated_vs_treated_upregulated_DEG_ID]
Shared_downregulated_DEGs_IDs<- pretreated_vs_treated_downregulated_DEG_ID[pretreated_vs_treated_downregulated_DEG_ID %in% untreated_vs_treated_downregulated_DEG_ID]

# Writing tables - Shared DEG IDs (Can't do p-values as they will be different for both calculations so Shared ID is the best we can do...)

write.table(Shared_DEGs_IDs, file = "DEGs/Shared/Shared_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_upregulated_DEGs_IDs, file = "DEGs/Shared/Shared_upregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)
write.table(Shared_downregulated_DEGs_IDs, file = "DEGs/Shared/Shared_downregulated_DEGs_IDs.txt", quote = F, row.names = F, col.names = F)

################## All counts ##################
