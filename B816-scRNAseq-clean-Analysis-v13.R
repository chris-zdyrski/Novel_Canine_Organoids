

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      ########   ##         ########        ##       ###     ##                                  
##     ##          ##         ##             ####      ####    ##                             
##    ##           ##         ##            ##  ##     ## ##   ##                                             
##    ##           ##         #######      ##    ##    ##  ##  ##                                  
##    ##           ##         ##          ##########   ##   ## ##                                   
##     ##          ##         ##          ##      ##   ##    ####                                   
##      ########   ########   ########    ##      ##   ##     ###                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ########      ########     #######    #######   ###     ##    ##########                                 
##     ##     ##     ##     ##     ##          ##          ##        ####    ##        ##                      
##     ##     ##     ##    ##      ##          ##          ##        ## ##   ##        ##                       
##     ########      #######       ########     #######    #######   ##  ##  ##        ##                               
##     ##            ##    ##      ##                 ##   ##        ##   ## ##        ##                       
##     ##            ##     ##     ##                 ##   ##        ##    ####        ##                       
##     ##            ##      ##    ########     #######    #######   ##     ###        ##                                
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
#FIGURES FOR PAPER (Seurat generic)


setwd("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2023-08 - UGA-RNAseq-COLLAB/2024-XX - Canine Organoid-RNAseq/05 scRNAseq - canine")

library("Seurat")
library("ggplot2")
library(dplyr)

#SAVE IN FORMS THAT LOAD FASTER
#load("FINAL_B816-organoids.RData")
#saveRDS(seurat_merged, file = "FINAL_B816-organoids.rds")
#qsave(seurat_merged, "FINAL_B816-organoids.qs")
seurat_merged <- qs::qread("FINAL_B816-organoids.qs")


all_cols<-scales::hue_pal()(13)
bladder_cols<-c("#F8766D","#E18A00","#BE9C00")
kidney_cols<-c("#8CAB00","#24B700","#00BE70")
lung_cols<-c("#00C1AB","#00BBDA","#8B93FF")
pancreas_cols<-c("#D575FE","#F962DD","#FF65AC")



umap_plot <- function(seurat_merged,color){
    DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
      ggtitle("B816 organoids: sub-types") + 
      theme(legend.position = "bottom") +
      scale_color_manual(values = color)
}

generate_feature_plot <- function(seurat_merged,gene="TOP2A") {
  #req(seurat_merged, gene)
  FeaturePlot(seurat_merged, features = gene) +
    scale_color_gradient(low = "grey", high = "red", name = paste0(gene, " (logCounts)")) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18)
    )
}

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#MERGED CELL FIGS:
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

umap_plot(seurat_merged,color=c(bladder_cols,kidney_cols,lung_cols,pancreas_cols))
generate_feature_plot (seurat_merged,gene="TOP2A")

DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
  ggtitle("B816 organoids: sub-types")


DimPlot(seurat_merged, reduction = "umap", group.by = "Tissue") +
  ggtitle("B816 organoids: sub-types")


tmp<-subset(seurat_merged, subset = Tissue == "Pancreas")
umap_plot(tmp,color=c(pancreas_cols))
generate_feature_plot (tmp,gene="TOP2A")

##################################################
#GENE-list provided by Chris
##################################################


DimPlot(tmp, reduction = "umap", group.by = "type2")+
  ggtitle("B816 organoids: pancreas")+
  scale_color_manual(values=c(pancreas_cols))

FeaturePlot(tmp, features = "TOP2A")


marker_genes_tf<-read.csv("05 CLEAN B816 - umap and Dotplot/scRNAseq genes.csv",as.is=T,header=T)

top_genes<-marker_genes_tf$Gene


dotplot <- DotPlot(seurat_merged, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6,group.by="type2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("DotPlot - All-Cells - Markers.pdf", plot = dotplot)


# (4) Heatmap of the top 10 genes per cell type


seurat_merged <- ScaleData(seurat_merged, features = top_genes)


heatmap <- DoHeatmap(seurat_merged, features = unique(top_genes), size = 3,group.by="type2") +
  ggtitle("Marker Genes per Cell Type")

ggsave("Heatmap - All-Cells - Markers.pdf", plot = heatmap)

####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#INDIVIDUAL SAMPLE SEURAT BOILER PLATE:
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


######################################
#MAKE INDIVIDUAL FILES
######################################

seurat_merged <- qs::qread("FINAL_B816-organoids.qs")
tmp_bladder <- subset(seurat_merged, subset = Tissue == "Bladder")
tmp_kidney <- subset(seurat_merged, subset = Tissue == "Kidney")
tmp_lung <- subset(seurat_merged, subset = Tissue == "Lung")
tmp_pancreas <- subset(seurat_merged, subset = Tissue == "Pancreas")


# Feature list for gene selection



######################################
#ALL FILES:
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_merged, reduction = "umap", group.by = "type2") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=all_cols)

seurat_obj<-seurat_merged

Idents(seurat_obj) <- seurat_obj@meta.data$type2

seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# Visualize the assigned cell types in the UMAP plot
DimPlot(seurat_obj, reduction = "umap", group.by = "type2") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=all_cols)


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("ALL_B816_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("All_B816_Heatmap.pdf", plot = heatmap)




######################################
#PANCREAS
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_pancreas, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid pancreatic cancer")+
  scale_color_manual(values=pancreas_cols)

seurat_obj<-tmp_pancreas

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Pancreas_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Pancreas_Heatmap.pdf", plot = heatmap)




######################################
#LUNG
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_lung, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid lung cancer")+coord_cartesian(xlim = c(-12, -4),ylim=c(-1.5,6))+
  scale_color_manual(values=lung_cols)

seurat_obj<-tmp_lung

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Lung_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Lung_Heatmap.pdf", plot = heatmap)










######################################
#BLADDER
######################################



# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_bladder, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid bladder cancer")+coord_cartesian(xlim = c(1, 8))+
  scale_color_manual(values=bladder_cols)

seurat_obj<-tmp_bladder

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Bladder_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Bladder_Heatmap.pdf", plot = heatmap)




######################################
#KIDNEY
######################################


# Visualize the assigned cell types in the UMAP plot
DimPlot(tmp_kidney, reduction = "umap", group.by = "cell_type") +
  ggtitle("Cell-Type:  organoid kidney cancer")+
  scale_color_manual(values=kidney_cols)

seurat_obj<-tmp_kidney

Idents(seurat_obj) <- seurat_obj@meta.data$cell_type


# (1) Generate a CSV file with the top genes per cell type
top_genes_per_cell_type <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)


# (2) Dotplot of the top few genes per cell type and their intensity/proportion
top_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  pull(gene)

dotplot <- DotPlot(seurat_obj, features = unique(top_genes), cols = c("blue", "red"), dot.scale = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top Genes per Cell Type")

ggsave("B816_Kidney_DotPlot.pdf", plot = dotplot)




# (4) Heatmap of the top 10 genes per cell type
heatmap_genes <- top_genes_per_cell_type %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  pull(gene)

seurat_obj <- ScaleData(seurat_obj, features = heatmap_genes)


heatmap <- DoHeatmap(seurat_obj, features = unique(heatmap_genes), size = 3) +
  ggtitle("Top 10 Genes per Cell Type")

ggsave("B816_Kidney_Heatmap.pdf", plot = heatmap)


