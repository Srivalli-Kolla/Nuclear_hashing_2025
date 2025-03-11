"""
Notebook to find Doublets and Demultiplex multinucleated data using HTODemux

**Created by :** Srivalli Kolla

**Created on :** 11 March, 2025

**Modified on :** 11 March, 2025

**University of Würzburg**

"""
  
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(cowplot)
library(tidyr)

#Load in data when you have following files: features.tsv.gz, barcodes.tsv.gz and matrix.mtx.gz
multiplexed.matrices = Read10X('../data/filtered_feature_bc_matrix/')

#view(multiplexed)
#str(multiplexed)

#Split the obtained list of two matrices in single matrices - one for **RNA expression** and one for **ADT expression**
multiplexed.RNA = multiplexed.matrices$`Gene Expression`
multiplexed.Antibodies = multiplexed.matrices$`Antibody Capture`

#Convert multiplexed.Antibodies into matrix
multiplexed.Antibodies = as.matrix(multiplexed.Antibodies)

include_list_hto = c('TotalSeqB1','TotalSeqB3','TotalSeqB4','TotalSeqB5','TotalSeqB6','TotalSeqB7','TotalSeqB8','TotalSeqB9')

multiplexed.HTO = subset(multiplexed.Antibodies, rownames(multiplexed.Antibodies) %in% include_list_hto)
length(include_list_hto)

#Create Seurat Object from RNA matrix and log-normalize the data
multiplexed <- CreateSeuratObject(counts = multiplexed.RNA)
multiplexed <- NormalizeData(multiplexed)

#Add the hashtag to each cell and normalize the data via CLR normalization
multiplexed[["HTO"]] <- CreateAssayObject(counts = multiplexed.HTO)
multiplexed <- NormalizeData(multiplexed, assay = "HTO", normalization.method = "CLR")

####Testing different quantiles for HTODemux()
multiplexed.test1 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.99) ## Default
FeatureScatter(multiplexed.test1, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3")
table(multiplexed.test1$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test1, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'../Demultiplexing_R/plots/ridgeplot_quant_0.99.png',base_height = 20, base_width = 20)      

multiplexed.test2 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.995)
FeatureScatter(multiplexed.test2, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3")
table(multiplexed.test2$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test2, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'../Demultiplexing_R/plots/ridgeplot_quant_0.995.png',base_height = 20, base_width = 20)

multiplexed.test3 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.95)
FeatureScatter(multiplexed.test3, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3")
table(multiplexed.test3$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test3, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'../Demultiplexing_R/plots/ridgeplot_quant_0.95.png',base_height = 20, base_width = 20)

multiplexed.test4 <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.9995)
FeatureScatter(multiplexed.test4, feature1 = "TotalSeqB1", feature2 = "TotalSeqB3")
table(multiplexed.test4$HTO_classification.global)
save_plot(plot = RidgePlot(multiplexed.test4, assay = "HTO", features = rownames(multiplexed[["HTO"]])),'../Demultiplexing_R/plots/ridgeplot_quant_0.9995.png',base_height = 20, base_width = 20)

#Observation from Ridge plots - I choose 0.9995 as it shows expression of particular hashtag nore than others and Alex took 0.995

#Demultiplex the samples
#Sometimes the positive quantile needs adjustment: Sometimes the positive.quantile needs to be adjusted: too stringent and you lose "good" cells. Not stringent enough and you don't remove enough doublets
multiplexed <- HTODemux(multiplexed, assay = "HTO", positive.quantile = 0.9995)

#Check the distribution of doublets/singlets/negative cells by looking at the Global classification results
table(multiplexed$HTO_classification.global)
save_plot(plot = barplot(table(multiplexed$HTO_classification.global)), '../Demultiplexing_R/plots/barplot_quant_0.9995.png')

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
multiplexed.negatives <-  RunTSNE(multiplexed.negatives, perplexity = 100)
save_plot(plot = DimPlot(multiplexed.negatives) ,'../Demultiplexing_R/plots/umap_singlet_doublets.png')

#Check for grouping of hashtags by tSNE
multiplexed.singlets <- subset(multiplexed,idents = 'Singlet' )
multiplexed.singlets <- FindVariableFeatures(multiplexed.singlets, selection.method = 'mean.var.plot')
multiplexed.singlets <-  ScaleData(multiplexed.singlets, features =  VariableFeatures(multiplexed.singlets))
multiplexed.singlets <-  RunPCA(multiplexed.singlets, features = VariableFeatures(multiplexed.singlets))
multiplexed.singlets <- FindNeighbors(multiplexed.singlets, reduction =  'pca', dims = 1:8)
multiplexed.singlets <-  FindClusters(multiplexed.singlets, resolution = 0.6, verbose =  FALSE)
multiplexed.singlets <- RunTSNE(multiplexed.singlets, reduction = 'pca', dims = 1:8)
save_plot(plot = DimPlot(multiplexed.singlets,group.by = "HTO_classification"), '../Demultiplexing_R/plots/umap_hashtags.png')

#Preparation for scRNAseq Analysis Workflow

#### Quality Control Filtering
#Keeping only the Singlets
multiplexed <- subset(multiplexed, idents = "Singlet")

#Add a column into the metadata for the content of mitochondrial transcripts
multiplexed[["percent.mt"]] = PercentageFeatureSet(multiplexed, pattern = "^mt-")

#Visualize QC metrics in violin plots-nFeature_RNA means xy Transcripts per Cell-nCount_RNA means xy Reads of Transcripts or total molecules per cell
VlnPlot(multiplexed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#First careful filtering: Remove cells with high mitochondrial transcript content (>10%) and high nCount_RNA
multiplexed <- subset(multiplexed, subset = percent.mt < 20 & nCount_RNA < 9000)

#Visualize as Violin Plot
VlnPlot(multiplexed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#How many cells are left after QC?
length(Cells(multiplexed))

#Indicate in the MetaData which cell comes from which sample

#Set Identities to HTO_classification
Idents(multiplexed) <- "HTO_classification"
levels(multiplexed)

#we can see that this is not the right order, so the first step is to fix the order in the Seurat object
levels(multiplexed) = c("TotalSeqB1","TotalSeqB3","TotalSeqB4","TotalSeqB5","TotalSeqB6" , "TotalSeqB7" ,"TotalSeqB8","TotalSeqB9")
levels(multiplexed)

#Assign sample names to the hashtags
sample.IDs = c("20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_10",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_11",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_12",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_13",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_14",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_15",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_16",
 "20241028_AC_Pkp2-MCMV_Mouse_Multiplex_Heart_17")

names(sample.IDs) <- levels(multiplexed)
multiplexed <- RenameIdents(multiplexed, sample.IDs)
multiplexed[["Sample"]] <- Idents(object = multiplexed)
levels(multiplexed)


#Indicate in the MetaData which cell comes from which experimental group

#Assign experimental condition names to the hashtags
exp.group.IDs = c("MCMV","MCMV","MCMV","MCMV","noninf","noninf","noninf","noninf")
names(exp.group.IDs) <- levels(multiplexed)
multiplexed <- RenameIdents(multiplexed, exp.group.IDs)
multiplexed[["Condition"]] <- Idents(object = multiplexed)
levels(multiplexed)

#Set Identities to Sample
Idents(multiplexed) <- "Sample"
levels(multiplexed)

A1_multiplexed.labels <- multiplexed@meta.data %>% 
  select(Sample, Condition)

# Saving files
write.csv(A1_multiplexed.labels, '../Demultiplexing_R/hashtags_demultiplex_R.csv')
save(multiplexed,file = "../Demultiplexing_R/hashtags_demultiplexing.Robj")
