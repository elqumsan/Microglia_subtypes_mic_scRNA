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
newObject$seurat_clusters <- NULL       # Remove this column to prevenr confusion 
newObject$cluster <- newObject[[paste0("RNA_snn_res.", optRes)]]

reorderCluster = c( "4", "9", "12", "13", "20", "14", "16",   # Prog
                    "17", "10", "18", "23",                   # Ery
                    "11", "24", "25", "7", "21", "28",         # Myeloid
                    "22","15","8","27", "26",              # B-Lymophoid
                    "5", "1", "0", "19", "6", "3", "2"     # T-Lymphoid
                    
                    )

newObject$cluster<- factor(newObject$cluster, levels = reorderCluster)
Idents(newObject) <- newObject$cluster # Set Seurat to use optRes

## Plot ideal resolution on tSNE and UMAP
nClust <- uniqueN(Idents(newObject))       # Setup color palette 
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired")) (nClust)
DimPlot(newObject, reduction = "umap", pt.size = 0.1 , label = TRUE, label.size = 3, cols = colCls ) + plotTheme + coord_fixed()


## Proportion / cell number composition per cluster
ggData <- as.matrix(prop.table(table(newObject$cluster, newObject$nCount_RNA), margin = 2))

pheatmap(ggData, color = colorRampPalette(colGEX)(10),
         cluster_rows = FALSE, cutree_cols = 2,
         display_numbers = TRUE, number_format = "%.3f", angle_col = 315,
         width = 4, height = 6, filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustComLibH.png" )

ggData <- data.frame(prop.table(table(newObject@meta.data$nCount_RNA, newObject$cluster),margin = 2))



colnames(ggData)<-c("library", "cluster", "value")

p1 <-  ggplot(ggData, aes(cluster, value, fill = library )) +
  geom_col() + xlab("cluster") + ylab("Proportion of cells (%)") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip() 
