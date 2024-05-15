################################################################################
### Load Libraries ###
################################################################################

library(Seurat)
library(SeuratObject)

### Create Seurat Object

VP4_SubLib1 <- c('VP4_SubLib1_Andinos', 'VP4_SubLib1_RD_ATCC')
dataset_loc <- "/Volumes/LVD_QVE-1/Projects/vp4dms/analysis/scRNAseq_analysis/Cell_Ranger/"

d10x.data <- sapply(VP4_SubLib1, function(i){
  d10x <- Read10X(file.path(dataset_loc,i,"outs/filtered_feature_bc_matrix"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x})
experiment.data <- do.call("cbind", d10x.data)

VP4.aggregate <- CreateSeuratObject(
  counts = experiment.data,
  project = "VP4",
  min.cells = 3,
  min.features = 10,
  names.field = 2,
  names.delim = "\\-")


head(VP4.aggregate[[]])

table(VP4.aggregate$orig.ident)
RidgePlot(VP4.aggregate, features="nFeature_RNA")


################################################################################
### Add % mito gene expression to seurat object metadata
################################################################################

VP4.aggregate[["percent.mt"]] <- PercentageFeatureSet(VP4.aggregate, pattern = "^MT-")
VP4.aggregate[["percent.ribo"]] <- PercentageFeatureSet(VP4.aggregate, pattern = "^RP[SL]")
VP4.aggregate[["percent.EVA71"]] <- PercentageFeatureSet(VP4.aggregate, pattern = "EVA71-Tainan")

VlnPlot(VP4.aggregate, features = "percent.mt", ncol = 1,  pt.size = 0.25) + 
  NoLegend() + ggtitle('Percent Mitochondrial Reads')

VlnPlot(VP4.aggregate, features = 'percent.ribo', ncol = 1, pt.size = 0.25) + 
  NoLegend() + ggtitle('Percent Ribosomal Reads')
VlnPlot(VP4.aggregate, features = 'percent.EVA71', ncol = 1, pt.size = 0.25) + 
  NoLegend() + ggtitle('Percent Viral Reads')

FeatureScatter(VP4.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  NoLegend() + facet_wrap(~VP4.aggregate$orig.ident)

plot4 <- FeatureScatter(VP4.aggregate, feature1 = "percent.mt", feature2 = "percent.ribo") + 
  NoLegend() + ggtitle('Percent Mitochondrial Reads vs. Percent Ribosomal Reads')+ 
  scale_y_log10() + scale_x_log10() + facet_wrap(~VP4.aggregate$orig.ident)
plot4

splitVP4 <- SplitObject(VP4.aggregate, split.by = "orig.ident")

################################################################################
### Perform integration analysis ###
################################################################################

for (i in VP4_SubLib1){
  
  halfmin_virus = min(splitVP4[[i]]$percent.EVA71[splitVP4[[i]]$percent.EVA71>0])
  
  percent.EVA71_log10 = log10(splitVP4[[i]]$percent.EVA71 + halfmin_virus)
  
  viruspercentage = mixtools::normalmixEM(splitVP4[[i]]$percent.EVA71[percent.EVA71_log10 > -2.5],k = 2)
  #call infected any cell above the mu + 2x sigma =  ----
  
  splitVP4[[i]]$InfectedStatus <- "Not_Infected"
  splitVP4[[i]]$InfectedStatus[splitVP4[[i]]$percent.EVA71 > (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)] <- "Infected"
  
  table(splitVP4[[i]]$InfectedStatus)
  
  splitVP4[[i]]$InfectedStatus_groups <- "Not_Infected"
  splitVP4[[i]]$InfectedStatus_groups[splitVP4[[i]]$percent.EVA71 > viruspercentage$mu[2]] <- "High"
  splitVP4[[i]]$InfectedStatus_groups[splitVP4[[i]]$percent.EVA71 < viruspercentage$mu[2]] <- "Low"
  splitVP4[[i]]$InfectedStatus_groups[splitVP4[[i]]$percent.EVA71 < (viruspercentage$mu[1]+viruspercentage$sigma[1]*2)]<- "Not_Infected" #annotate the threshold in the meta.data object
}

plot_list = list()
for (i in VP4_SubLib1){
  plot_list[[i]] = VlnPlot(splitVP4[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1, group.by = "InfectedStatus") +ggtitle(i) +theme(plot.title = element_text(hjust = 1, vjust=4.2))
}
plot_list

VP4_SubLib1_RD_Andino <- splitVP4[['VP4_SubLib1_Andinos']]

VP4_SubLib1_RD_ATCC <- splitVP4[['VP4_SubLib1_RD_ATCC']]

doublets <- read.table("/Volumes/lvd_qve/Projects/vp4dms/analysis/scRNAseq_analysis/Seurat/VP4_SubLib1_RD_ATCC/VP4_SubLib1_RD_ATCC_Doublet_scores.tsv", header = F)
colnames(doublets) <- c("doublet_scores","predicted_doublets")

VP4_SubLib1_RD_ATCC <- AddMetaData(VP4_SubLib1_RD_ATCC, metadata = doublets$doublet_scores, col.name = "doublet_scores")
VP4_SubLib1_RD_ATCC <- AddMetaData(VP4_SubLib1_RD_ATCC, metadata = doublets$predicted_doublets, col.name = "predicted_doublets")
head(VP4_SubLib1_RD_ATCC[[]])
table(VP4_SubLib1_RD_ATCC$predicted_doublets)

doublets <- read.table("/Volumes/LVD_QVE-1/Projects/vp4dms/analysis/scRNAseq_analysis/Cell_Ranger/VP4_SubLib1_Andinos/VP4_SubLib1_Andinos_Doublet_scores.tsv", header = F)
colnames(doublets) <- c("doublet_scores","predicted_doublets")

VP4_SubLib1_RD_Andino <- AddMetaData(VP4_SubLib1_RD_Andino, metadata = doublets$doublet_scores, col.name = "doublet_scores")
VP4_SubLib1_RD_Andino <- AddMetaData(VP4_SubLib1_RD_Andino, metadata = doublets$predicted_doublets, col.name = "predicted_doublets")
head(VP4_SubLib1_RD_Andino[[]])
table(VP4_SubLib1_RD_Andino$predicted_doublets)

VP4_SubLib1_RD_Andino <-subset(VP4_SubLib1_RD_Andino, subset = nCount_RNA > 10000 & nFeature_RNA >4000 & percent.mt < 20 & doublet_scores < 0.608)

VP4_SubLib1_RD_ATCC <- subset(VP4_SubLib1_RD_ATCC, subset = nCount_RNA > 10000 & nFeature_RNA >4000 & percent.mt < 20 & doublet_scores <0.489)

table(VP4_SubLib1_RD_Andino$predicted_doublets)
table(VP4_SubLib1_RD_ATCC$predicted_doublets)
###Subset by sample according to VlnPlots and doublet scores

seurat_object_merged_filtered_list <- merge(VP4_SubLib1_RD_ATCC, y = VP4_SubLib1_RD_Andino)

seurat_object_merged_filtered_list <- SplitObject(seurat_object_merged_filtered_list, split.by = 'orig.ident')
# normalize each dataset individually and find 2000 variable features
seurat_object_merged_filtered_list <- lapply(X = seurat_object_merged_filtered_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_object_merged_filtered_list)
seurat_object_merged_filtered_list <- lapply(X = seurat_object_merged_filtered_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = seurat_object_merged_filtered_list,
                                  anchor.features = features,
                                  reduction = "rpca",
                                  dims = 1:50)


# Integrate datasets
seurat_object_integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

################################################################################
### Cell cycle scoring, dimensionality reduction and clustering ###
################################################################################
exp.mat <- read.table(file = "/Users/mariskanishca/Downloads/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE,
                      as.is = TRUE, row.names = 1)

seurat_object_integrated <- CellCycleScoring(seurat_object_integrated, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

DefaultAssay(seurat_object_integrated) <- "integrated"

seurat_object_integrated <- ScaleData(seurat_object_integrated, verbose = FALSE)
seurat_object_integrated <- RunPCA(seurat_object_integrated, npcs = 30, verbose = FALSE)
seurat_object_integrated <- RunUMAP(seurat_object_integrated, reduction = "pca", dims = 1:30)
seurat_object_integrated <- FindNeighbors(seurat_object_integrated, reduction = "pca", dims = 1:30)

# cluster cells from 0.1 to 1 resolution for later exploration
for(i in 1:length(seq(0.1, 1, 0.1))){
  
  seurat_object_integrated <- FindClusters(seurat_object_integrated, resolution = seq(0.1, 1, 0.05)[i])
  
}
seurat_object_integrated <- FindClusters(seurat_object_integrated, resolution = 0.15)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.1', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.15', label =TRUE)
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.2', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.25', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.3', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.35', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.4', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.45', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.5', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.55', label =TRUE, split.by = 'orig.ident')
DimPlot(seurat_object_integrated, group.by = 'integrated_snn_res.0.55', label =TRUE, split.by = 'orig.ident')

FeaturePlot(seurat_object_integrated, features = c('percent.EVA71', 'percent.mt'), split.by = 'orig.ident')
clustree(seurat_object_integrated, assay = 'integrated')
################################################################################
### Dimensionality reduction and clustering of RNA assay ###
################################################################################
Idents(seurat_object_integrated) <- 'integrated_snn_res.0.15'
DefaultAssay(seurat_object_integrated) <- 'RNA'

FeaturePlot(seurat_object_integrated, features = c("OASL","IFIT2","ISG15","CCL5","IFIT1","ISG20"), split.by = "orig.ident")


seurat_object_integrated@meta.data$condition<-paste0(seurat_object_integrated@meta.data$orig.ident)
Idents(seurat_object_integrated)<-"condition"

df_contrasts<-t(as.data.frame(combn(levels(as.factor(seurat_object_integrated@meta.data[["condition"]])),2)))
colnames(df_contrasts)<-c("d1","d2")
rownames(df_contrasts)<-c(1:nrow(df_contrasts))
df_contrasts<-as.data.frame(df_contrasts)
##############
#####select the relevant rows for your analysis here#######
df_contrasts_selec<-df_contrasts[c(1),]
##############
##############
DE_between_condition<-vector(mode = "list", length = nrow(df_contrasts_selec))
DE_between_condition_upregulated<-vector(mode = "list", length = nrow(df_contrasts_selec))
DE_between_condition_downregulated<-vector(mode = "list", length = nrow(df_contrasts_selec))
#############
#############
for (x in c(1:nrow(df_contrasts_selec))) {
  print(x)
  
  cell.count1<-colnames(seurat_object_integrated)[which(seurat_object_integrated@meta.data$condition==df_contrasts_selec$d1[x])]
  
  print(length(cell.count1))
  
  cell.count2<-colnames(seurat_object_integrated)[which(seurat_object_integrated@meta.data$condition==df_contrasts_selec$d2[x])]
  
  print(length(cell.count2))
  
  min.cell.count<-min(length(cell.count1),length(cell.count2))
  
  if (min.cell.count>=50) {
    
    DE_between_condition[[x]]<-FindMarkers(object = seurat_object_integrated, ident.1 = df_contrasts_selec$d1[x], ident.2= df_contrasts_selec$d2[x], min.pct = 0.10, logfc.threshold = 0.01, min.cells.group=50,test.use = "wilcox")
    
    DE_between_condition[[x]]$Gene<-rownames(DE_between_condition[[x]])
    
    DE_between_condition[[x]]$FC<-2^(DE_between_condition[[x]]$avg_log2FC)
    
    DE_between_condition_upregulated[[x]]<-DE_between_condition[[x]][which(DE_between_condition[[x]]$FC>1),]
    DE_between_condition_downregulated[[x]]<-DE_between_condition[[x]][which(DE_between_condition[[x]]$FC<1),]
  }
  else {
    DE_between_condition_upregulated[[x]]<-data.frame(Result=c("not_enough_cells_to_compare"))
    DE_between_condition_downregulated[[x]]<-data.frame(Result=c("not_enough_cells_to_compare"))
  }
  
}

##############
##############

###Data visualization for DE - volcano plot

ggplot(DE_between_condition[[1]])+
  geom_point(aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
  ggrepel::geom_text_repel(aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene))

################
################
library(ggrepel)
joined_seurat.markers$delabel <- NA
joined_seurat.markers$diffexpressed <- "NO"
# if log2Foldchange > 2 and pvalue < 0.05, set as "UP" 

joined_seurat.markers$diffexpressed[joined_seurat.markers$avg_log2FC > 2 & joined_seurat.markers$p_val_adj < 0.05] <- "UP"

# if log2Foldchange < -2 and pvalue < 0.05, set as "DOWN"
joined_seurat.markers$diffexpressed[joined_seurat.markers$avg_log2FC < -2 & joined_seurat.markers$p_val_adj < 0.05] <- "DOWN"

joined_seurat.markers$delabel[joined_seurat.markers$diffexpressed != "NO"] <- joined_seurat.markers$gene[joined_seurat.markers$diffexpressed != "NO"]

#Create a column for the shape manual
joined_seurat.markers$clusters <- factor(joined_seurat.markers$cluster)

p1 <- ggplot(joined_seurat.markers, aes(avg_log2FC, -log(p_val_adj,10), shape = factor(cluster))) + # -log10 conversion 
  geom_point(size=1) +
  xlab(expression("log"[2]*"(Fold Change)")) + 
  ylab(expression("-log"[10]*"Pvalue"))+
  geom_vline(xintercept=c(-2, 2), col="grey", linetype = 'dashed') +
  geom_hline(yintercept=-log10(0.05), col="grey", linetype = 'dashed')+
  coord_cartesian(ylim = c(0, 300), xlim = c(-5,5))+
  geom_point(data=subset(joined_seurat.markers, diffexpressed == "UP"), color="darkred")+
  geom_point(data=subset(joined_seurat.markers, diffexpressed == "DOWN"), color="darkblue")+
  #scale_color_brewer(palette="Paired") + theme(
  #legend.position = c(1, 0),
  #legend.justification = c("right", "bottom"),
  #legend.box.just = "right",
  #legend.margin = margin(6, 6, 6, 6)
  #)+
  labs(shape = 'Clusters') + geom_point(data=subset(joined_seurat.markers, diffexpressed == "NO"), color="grey")+
  scale_shape_manual(values=1:nlevels(joined_seurat.markers$clusters))+
  facet_wrap(~factor(cluster)) + labs(title = "Differential Expression across Clusters") +NoLegend()

p1 + geom_label_repel(size = 2, aes(label = delabel), max.overlaps = 20)

##############
#Use EnhancedVolcano to produce a volcano plot
##############
EnhancedVolcano(DE_between_condition[[1]], lab = rownames(DE_between_condition[[1]]), x= 'avg_log2FC', y= 'p_val_adj',  title = 'ATCC RD Cells versus Andino RD Cells', FCcutoff = 0.5, xlim = c(-4,4), subtitle = NULL, caption = NULL)
