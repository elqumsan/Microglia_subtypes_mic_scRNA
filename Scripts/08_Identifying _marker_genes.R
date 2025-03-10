library(clustree)
library("RColorBrewer")
library(ggplot2)
library(gridExtra)

library(data.table)
library(Matrix)
library(ggplot2)
library(plotly)
library(patchwork)
library(RColorBrewer)
library(Seurat)
library(reticulate)
library(pdist)
library(phateR)
library(SeuratWrappers)
library(monocle3)
library(slingshot)
library(tradeSeq)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(BiocParallel)

#### Define color palettes and plot themes
plotTheme <- theme_classic(base_size = 18)
colLib = brewer.pal(7,  "Paired")
names(colLib)= c("BM0pos", "BM0neg", "BM7pos", "BM7eng", "BM9pos", "BM9neg", "BM2pos")
colDnr = colLib[c(1,3,5,7)]
names(colDnr) = c("BM0", "BM7", "BM9", "BM2")
colGEX <- c("grey85", brewer.pal(7, "Reds"))
colCcy = c("black", "blue", "darkorange")
os <- import("os")
##### Perform Clustering

newObject <- FindNeighbors(integrated.strain, dims = 1:20)
newObject <- FindClusters(newObject, resolution = seq(0.4, 1.5, 0.1))

#### Find ideal resolution
 p1 <-clustree::clustree(newObject)

 

ggsave(p1 + scale_color_manual(values = brewer.pal(12, "Paired")) +
         guides(color= guide_legend(ncol = 2)),
       width = 10, height = 8, filename= "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustClustree.png"
       )

optRes <- 1.0                           # determined from clustering tree
newObject$seurat_clusters <- NULL       # Remove this column to prevent confusion 
newObject$cluster <- newObject[[paste0("RNA_snn_res.", optRes)]]

reorderCluster = c( "4", "9", "12", "13", "20", "14", "16",   # Prog
                    "17", "10", "18", "23",                   # Ery
                    "11", "24", "25", "7", "21",             # Myeloid
                    "22","15","8", "26",                     # B-Lymophoid
                    "5", "1", "0", "19", "6", "3", "2"     # T-Lymphoid
                    
                    )

newObject$cluster<- factor(newObject$cluster, levels = reorderCluster)
Idents(newObject) <- newObject$cluster # Set Seurat to use optRes

## Plot ideal resolution on tSNE and UMAP
nClust <- uniqueN(Idents(newObject))       # Setup color palette 
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired")) (nClust)
p1 <- DimPlot(newObject, reduction = "tsne", pt.size = 0.1,  label = TRUE, label.size = 3, cols = colCls) + plotTheme + coord_fixed()

p2 <- DimPlot(newObject, reduction = "umap", pt.size = 0.1 , label = TRUE, label.size = 3, cols = colCls ) + plotTheme + coord_fixed()

ggsave(p1+ p2+ plot_layout(guides = "collect"), 
       width = 10, height = 4, filename  = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clusDimPlot.png")

## Proportion / cell number composition per sub_type
ggData <- as.matrix(prop.table(table(newObject$cluster, newObject$orig.ident), margin = 2))

pheatmap(ggData, color = colorRampPalette(colGEX)(100),
         cluster_rows = FALSE, cutree_cols = 2,
         display_numbers = TRUE, number_format = "%.3f", angle_col = 315,
         width = 4, height = 6, filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustComLibH.png" )



ggData <- data.frame(prop.table(table(newObject$cluster, newObject$orig.ident), margin = 2))

colnames(ggData)<-c("cluster.No", "sub_type", "proportion")

p1 <-  ggplot(ggData, aes(sub_type, proportion, fill = cluster.No )) +
  geom_col() + xlab("sub_type") + ylab("Proportion of cells (%)") +
  scale_fill_manual(values = colCls) + plotTheme + coord_flip() 

ggData <- data.frame(table(newObject$cluster, newObject$orig.ident))
colnames(ggData)<- c("cluster.No", "sub_type", "proportion")

p2 <- ggplot(ggData, aes(sub_type, proportion, fill = cluster.No))+
       geom_col() + xlab("sub_type") + ylab("Cell Number") + 
        scale_fill_manual(values = colCls) + plotTheme + coord_flip()


ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 14, height = 6 , filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustComLib.png" )


####### Proportion / cell number composition per cluster
ggData= data.frame(prop.table(table(newObject$orig.ident, newObject$cluster), margin = 2))
colnames(ggData) <- c("sub_type", "cluster.No", "proportion")
p1 <- ggplot(ggData, aes(cluster.No, proportion , fill= sub_type)) + 
     geom_col() + xlab("cluster.No") + ylab("proportion of cells (%)") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip()

ggData = data.frame(table(newObject$orig.ident, newObject$cluster))
colnames(ggData) = c("sub_type", "cluster.No", "proportion")
p2 <- ggplot(ggData, aes(cluster.No, proportion, fill = sub_type)) +
  geom_col() + xlab("cluster") + ylab(" Cell Number") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip()

ggsave(p1 + p2 + plot_layout(guides = "collect"),
       width = 10, height = 6,  
       filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustComClust.png")


#########################
######################### Find Markers

oupMarkers <- FindAllMarkers(newObject, only.pos = TRUE , logfc.threshold =1.0,  min.pct = 0.2)
oupMarkers <- data.table(oupMarkers)
oupMarkers$pct.diff = oupMarkers$pct.1 = oupMarkers$pct.2 
oupMarkers <-oupMarkers[,c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]
fwrite(oupMarkers, sep = "\t", file = "../Microglia_subtypes_mic_scRNA/findings/06_differential_expression_analysis/clusterMarkers.txt")


#### Check if known genes are  in the marker gene list
