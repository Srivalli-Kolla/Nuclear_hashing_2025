# Notebook to find Doublets and Demultiplex multinucleated data using HTODemux
# 
# Created by : Srivalli Kolla
# Created on : 11 March, 2025
# Modified on : 01 April, 2025
#
# University of Würzburg

library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(cowplot)
library(tidyr)

#Load in data when you have h5ad
multiplexed.matrices = zellkonverter::readH5AD('Documents/GitHub/Nuclear_hashing_2025/data/0.01_full/cb_0.01_full_with_hashtags_06_03_25.h5ad')

# Convert to Seurat object
seurat_obj <- as.Seurat(multiplexed.matrices, counts = "X", data = "X")

# Extract RNA and Antibody Capture (Hashtag) data
multiplexed.Antibodies <- GetAssayData(seurat_obj)

#Convert multiplexed.Antibodies into matrix
multiplexed.Antibodies = as.matrix(multiplexed.Antibodies)

include_list_hto = c('TotalSeqB1','TotalSeqB3','TotalSeqB4','TotalSeqB5','TotalSeqB6','TotalSeqB7','TotalSeqB8','TotalSeqB9')

multiplexed.HTO = subset(multiplexed.Antibodies, rownames(multiplexed.Antibodies) %in% include_list_hto)
length(include_list_hto)

#Create Seurat Object from RNA matrix and log-normalize the data
multiplexed <- CreateSeuratObject(counts = multiplexed.Antibodies)
multiplexed <- NormalizeData(multiplexed)

#Add the hashtag to each cell and normalize the data via CLR normalization
multiplexed[["HTO"]] <- CreateAssayObject(counts = multiplexed.HTO)
multiplexed <- NormalizeData(multiplexed, assay = "HTO", normalization.method = "CLR")

####Testing different quantiles for HTODemux()
multiplexed.test1 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.99) ## Default
save_plot(plot = FeatureScatter(multiplexed.test1, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3"),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/scatterplot_quant_0.99_cb_filtered.png')
table(multiplexed.test1$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test1, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/ridgeplot_quant_0.99_cb_filtered.png',base_height = 20, base_width = 20)      

multiplexed.test2 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.995)
save_plot(plot = FeatureScatter(multiplexed.test2, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3"),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/scatterplot_quant_0.995_cb_filtered.png')
table(multiplexed.test2$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test2, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/ridgeplot_quant_0.995_cb_filtered.png',base_height = 20, base_width = 20)

multiplexed.test3 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.95)
save_plot(plot = FeatureScatter(multiplexed.test3, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3"),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/scatterplot_quant_0.95_cb_filtered.png')
table(multiplexed.test3$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test3, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/ridgeplot_quant_0.95_cb_filtered.png',base_height = 20, base_width = 20)

multiplexed.test4 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.9995)
save_plot(plot = FeatureScatter(multiplexed.test4, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3"),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/scatterplot_quant_0.9995_cb_filtered.png')
table(multiplexed.test4$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test4, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/ridgeplot_quant_0.9995_cb_filtered.png',base_height = 20, base_width = 20)

#Observation from Ridge plots - I choose 0.95 as it shows expression of particular hashtag nore than others and Alex took 0.995

#Demultiplex the samples
#Sometimes the positive quantile needs adjustment: Sometimes the positive.quantile needs to be adjusted: too stringent and you lose "good" cells. Not stringent enough and you don't remove enough doublets
multiplexed <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.95)

#Check the distribution of doublets/singlets/negative cells by looking at the Global classification results
table <- table(multiplexed$HTO_classification.global)
table
png("Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/barplot_quant_0.9995_cb_filtered.png", width = 800, height = 600)
barplot_heights <- barplot(table, ylab = "Count", las = 2, main = "HTO Classification")
text(x = barplot_heights, y = table, labels = table, pos = 1)
dev.off()

#Check the number of unique RNA molecules in each cell according to each hashtag status
Idents(multiplexed) = "HTO_maxID"
RidgePlot(multiplexed, assay = "HTO", features = rownames(multiplexed[["HTO"]]))
VlnPlot(multiplexed, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

#FeatureScatter plots to see if positive.quantile is well chosen
FeatureScatter(multiplexed, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3")
FeatureScatter(multiplexed, feature1 = "TotalSeqB4", feature2 = "TotalSeqB5")
FeatureScatter(multiplexed, feature1 = "TotalSeqB6", feature2 = "TotalSeqB7")
FeatureScatter(multiplexed, feature1 = "TotalSeqB8", feature2 = "TotalSeqB9")

#Check the global Hashtag signal in all cells
HTOHeatmap(multiplexed, assay = "HTO", ncells = 10000)

#Check for grouping of singlets and doublets by calculating tSNE
Idents(multiplexed) <- "HTO_classification.global" ## Toi know singlet, doublet, negative status
multiplexed.negatives <- subset(multiplexed,idents = 'Negative', invert = TRUE)
DefaultAssay(multiplexed.negatives) <-  'HTO'
multiplexed.negatives <- ScaleData(multiplexed.negatives, features = rownames(multiplexed.negatives))
multiplexed.negatives <-  RunPCA(multiplexed.negatives, features = rownames(multiplexed.negatives),approx = FALSE)
#multiplexed.negatives <-  RunTSNE(multiplexed.negatives, perplexity = 100)
save_plot(plot = DimPlot(multiplexed.negatives) ,'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/umap_singlet_doublets_cb_filtered.png')

#Check for grouping of hashtags by tSNE
Idents(multiplexed) <- "HTO_classification.global"
multiplexed.clustering <- subset(multiplexed,idents = 'Singlet')
multiplexed.clustering <-  FindVariableFeatures(multiplexed.clustering,selection.method = 'mean.var.plot')
multiplexed.clustering <-  ScaleData(multiplexed.clustering, features = VariableFeatures(multiplexed.clustering))
multiplexed.clustering <- RunPCA(multiplexed.clustering)
multiplexed.clustering <-  RunTSNE(multiplexed.clustering, dims = 1:9, perplexity = 100)
save_plot(plot = DimPlot(multiplexed.clustering,group.by = "HTO_classification"), 'Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/umap_hashtags_cb_filtered.png')

#Check the distribution of doublets/singlets/negative cells by looking at the Global classification results
table <- table(multiplexed.clustering$HTO_classification)
table
png("Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/plots/barplot_quant_0.9995_hash_cb_filtered.png")
par(mar = c(8, 4, 4, 2)) 
barplot_heights <- barplot(table, ylab = "Count", main = "HTO Classification", names.arg = names(table), las = 2)
text(x = barplot_heights, y = table, labels = table, pos = 3, cex = 0.8)
dev.off()

# Convert Seurat to SingleCellExperiment to save the file
sce <- as.SingleCellExperiment(multiplexed)
writeH5AD(sce, "Documents/GitHub/Nuclear_hashing_2025/Demultiplexing_R/demultiplexed_HTODemux_cb_filtered.h5ad")