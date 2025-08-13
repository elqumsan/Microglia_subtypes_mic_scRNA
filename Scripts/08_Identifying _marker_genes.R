# install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")

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
library(slingshot)
library(SingleCellExperiment)
library(BiocParallel)
library(mgcv)

library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(BiocParallel)
library(msigdb)
library(msigdbr)
library(ExperimentHub)
library(GSEABase)
library(enrichR)

library(monocle3)
library(tradeSeq)

library(gProfileR)

library(org.Mm.eg.db)
library(DOSE)

#### Define color palettes and plot themes
plotTheme <- theme_classic(base_size = 18)
colLib = brewer.pal(2,  "Paired")
names(colLib)= c("Veh", "AZT")
colDnr = colLib[c(1,3,5,7)]
names(colDnr) = c("BM0", "BM7", "BM9", "BM2")
colGEX <- c("grey85", brewer.pal(7, "Reds"))
colCcy = c("black", "blue", "darkorange")


colorMedia = c("black", "blue", "magenta")
names(colorMedia) = c("D21-fm", "D21-nr" , "D21-tr")

os <- import("os")
##### Perform Clustering

newObject <- FindNeighbors(integrated.strain, dims = 0:12)
newObject <- FindClusters(newObject, resolution = seq( 1.0 ,1.5, 0.1))

#### Find ideal resolution
 p1 <-clustree::clustree(newObject)

 

ggsave(p1 + scale_color_manual(values = brewer.pal(12, "Paired")) +
         guides(color= guide_legend(ncol = 2)),
       width = 10, height = 8, filename= "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustClustree.png"
       )

optRes <- 1.0                           # determined from clustering tree
newObject$seurat_clusters <- NULL       # Remove this column to prevent confusion 
newObject$cluster <- newObject[[paste0("RNA_snn_res.", optRes)]]


#reorderCluster = c( "4", "9", "12", "13", "20", "14", "16",   # Prog
#                    "17", "10", "18", "23",                   # Ery
#                    "11", "24", "25", "7", "21",             # Myeloid
#                    "22","15","8", "26",                     # B-Lymophoid
#                    "5", "1", "0", "19", "6", "3", "2"     # T-Lymphoid
#                    
#                   )

reorderCluster = c( "0", "1", "2","3","4","5", "6", "7", "8", "9", "10", "11", "12" )

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
  scale_fill_manual(values = brewer.pal(9,  "Paired")) + theme_classic(base_size = 18) +  coord_flip()

ggData = data.frame(table(newObject$orig.ident, newObject$cluster))
colnames(ggData) = c("sub_type", "cluster.No", "proportion")
p2 <- ggplot(ggData, aes(cluster.No, proportion, fill = sub_type)) + 
  geom_col() + xlab("cluster") + ylab(" Cell Number") +
  scale_fill_manual(values = brewer.pal(9,  "Paired")) + theme_classic(base_size = 18) + coord_flip() 

ggsave(p1 + p2 + plot_layout(guides = "collect"),
       width = 10, height = 6,  
       filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/clustComClust.png")


#########################
######################### Find Markers genes

oupMarker <- FindAllMarkers(newObject, only.pos = TRUE , logfc.threshold =1.0,  min.pct = 0.2)
oupMarker <- data.table(oupMarker)
oupMarker$pct.diff = oupMarker$pct.1 - oupMarker$pct.2 
oupMarker <-oupMarker[,c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "pct.diff", "p_val", "p_val_adj")]


fwrite(oupMarker, sep = "\t", file = "../Microglia_subtypes_mic_scRNA/findings/06_differential_expression_analysis/clusterMarkers.txt")


#### Check if known genes are  in the marker gene list
KnownGenes <- c("Tmem119", "P2ry12", "Cx3cr1", "Trem2", "Sall1", "Aif1","Iba1", "Apoe", "Lpl", "Trem2", "Cst7" , "Itgax","CD11c", "Ifitm3",
                "B2m","Siglech", "Pu.1", "Irf8", "Sall1", "Runx1")
KnownGenes <- c("Cst7",  "Apoe", "Cx3cr1", "Tmem119", "Ifit3", "Ifitm3", "Irf7", "Hexb", "Cd81", "Cst3", "Rplp1", "Rps21", "Rps24" 
                , "C3ar1", "Stmn1", "Top2a", "Birc5")
oupMarker[gene %in% KnownGenes]

## Get top genes for each cluster and do dotplot / violin plot
oupMarker$cluster = factor(oupMarker$cluster, levels = reorderCluster)
oupMarker = oupMarker[order(cluster, -avg_log2FC)]
genes.to.plot <- unique(oupMarker[cluster %in% reorderCluster,
                                  head(.SD, 2), by = "cluster"]$gene)

p1 <- DotPlot(newObject, group.by = "cluster", features = genes.to.plot) + 
  coord_flip() + scale_color_gradientn(colors = colGEX)  +
                                         theme(axis.text.x = element_text(angle = - 45, hjust = 0) )
ggsave(p1, width = 10, height = 8, filename = "../Microglia_subtypes_mic_scRNA/findings/06_differential_expression_analysis/clusterMarkersDot.png" )

p2 <- VlnPlot(newObject, group.by = "cluster", fill.by = "ident", cols= colCls, features = genes.to.plot, stack = TRUE, flip = TRUE )

ggsave(p2, width = 10, height = 8 , filename = "../Microglia_subtypes_mic_scRNA/findings/06_differential_expression_analysis/clustMarkersVln.png")

####### Gene module analysis + module score
####  Functional analysis using clusterProfiler and  msigdb

oupMarkerFunc = data.table()
set.seed(42)

for (iDB in c("C8", "C5_GO:BP")) {
# Get reference gene set from msigdbr  
  inpGS <- tstrsplit(iDB, "_")
  msigCat <- inpGS[[1]] ;  msigSubCat <- NULL
  if (length(inpGS) >= 2) { msigSubCat <- inpGS[[2]]}
  
#  inpGS <- data.frame(msigdbr( species = "mouse", category = "C2", subcategory = "CGP" ))
  inpGS <- data.frame(msigdbr( species = "Mus musculus", collection = msigCat, subcollection = msigSubCat  ))
  
  #  inpGS <- data.frame(msigdbr( species = "Mus musculus", category ="M2", subcategory = "CGP" ))
  
  
  inpGS <- inpGS[, c("gs_name", "gene_symbol")]
  
  #### start up clusterProfiler
  for (i in unique(oupMarker$cluster)) {
  tmpOut <-  enricher(oupMarker[cluster == i]$gene, universe = rownames(newObject) , TERM2GENE = inpGS)
    tmpOut <- data.table(sigdb = iDB, cluster = i , data.frame(tmpOut[, -2])) 
    tmpOut$mLog10Padj <- -log10(tmpOut$p.adjust)
    tmpOut <- tmpOut[order(-mLog10Padj)]
    oupMarkerFunc <- rbindlist(list(oupMarkerFunc, tmpOut))
      }
  
} ### end of iDB for loop


oupMarkerFunc$cluster <- factor(oupMarkerFunc$cluster, levels = reorderCluster)
newObject@misc$markerFunc <- oupMarkerFunc ## store func analysis into Seurat object 

#### plot functional analysis results
ggData <- oupMarkerFunc[grep("MANNO_MIDBRAIN_NEUROTYPES_HRGL1",ID)]
ggData$ID <- gsub("brain-resident macrophages", "", ggData$ID)
ggData$ID <- factor(ggData$ID, levels = unique(ggData$ID))

p1 <- ggplot(ggData, aes(cluster, ID, size = Count, color = mLog10Padj)) +
  geom_point() + ggtitle("Brain-resident macrophages")+
  theme_linedraw(base_size = 18) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  scale_color_gradientn(colours = colGEX, limits = c(0, 10), na.value = colGEX[8])

ggsave(p1, width = 12,height = 8, filename = "../Microglia_subtypes_mic_scRNA/findings/08_identifying_marker_genes/clustMarkersFuncMB.png")



tmp <- oupMarkerFunc[sigdb == "C5_GO:BP"]
tmp <-tmp[grep("DIFF", ID)][!grep("POSITIVE|NEGATIVE", ID)]$ID
ggData <- oupMarkerFunc[ID %in% tmp]
substr(gsub("GOBP_", "", ggData$ID), 1 , 30)
ggData$ID <- factor(ggData$ID, levels = unique(ggData$ID))

p1 <- ggplot(ggData, aes(cluster, ID, size= Count, color = mLog10Padj)) + 
  geom_point() + ggtitle(" DIFF signatures") + 
  theme_linedraw(base_size = 18) +
  theme(axis.text.x = element_text(angle =-45, hjust = 0 )) +
  scale_colour_gradientn(colours = colGEX, limits = c(0, 8), na.value = colGEX[8])

ggsave(p1, width = 12 , height = 10 , filename= "../Microglia_subtypes_mic_scRNA/findings/08_identifying_marker_genes/clustMarkersFuncDiff.png")


##### cell cycle + plot on UMAP
newObject <- CellCycleScoring(newObject, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes )
p1 <- DimPlot(newObject, reduction = "umap", pt.size = 0.1, group.by = "cluster",
              shuffle = TRUE , cols = colCls ) + plotTheme + coord_fixed()

p2 <- DimPlot(newObject, reduction = "umap", pt.size = 0.1 , group.by = "Phase", 
        shuffle = TRUE, cols = colCcy ) + plotTheme + coord_fixed()

p3 <- FeaturePlot(newObject, pt.size = 0.1,  feature ="S.Score") + 
      scale_color_distiller( palette = "RdY1Bu") + plotTheme + coord_fixed()

p4 <- FeaturePlot(newObject, pt.size = 0.1, feature = "G2M.Score") + scale_color_distiller( palette = "RdY1Bu") +
     plotTheme + coord_fixed()

ggsave(p1 + p2 + p3 + p4 , width = 10, height = 8 , 
       filename = "../Microglia_subtypes_mic_scRNA/findings/08_identifying_marker_genes/clusModuCellCycle.png")


### Add module score (Here I use using HallMark Signaling gene sets)
inpGS <- data.table(msigdbr(species = "Mus musculus" , collection = "C8"))
inpGS <- inpGS[grep("SIGNAL", gs_name)]
inpGS$gs_name <- gsub("HALLMARK", "", inpGS$gs_name)
inpGS$gs_name <- gsub("_SIGNALING" , "", inpGS$gs_name)
inpGS <- split(inpGS$gene_symbol , inpGS$gs_name)

newObject <- AddModuleScore(newObject, feature = inpGS , name = "HALLMARK")

colnames(newObject@meta.data) [grep("HALLMARK" , colnames(newObject@meta.data))] <- names(inpGS)

p1 <- FeaturePlot(newObject, reduction = "umap", pt.size = 0.1 , features = names(inpGS), order = TRUE ) &
  scale_color_distiller(palette = "RdY1Bu") & plotTheme & coord_fixed()

ggsave(p1, width = 20 , height = 12, filename= "../Microglia_subtypes_mic_scRNA/findings/08_identifying_marker_genes/clustModuScore.png")


###### Annotating clusters

oupAnnot = data.table(cluster = reorderCluster)
## Add in top 5 marker genes
ggData <- oupMarker[, head(.SD, 8), by = "cluster"]
ggData <- ggData[,paste0(gene, collapse = ","), by = "cluster"]
colnames(ggData)[2] <- "markers"
oupAnnot <- ggData[oupAnnot, on = "cluster"]

### Add in Microglia_Homeostatic enrichment
ggData <- oupMarkerFunc[grep("BONE_MARROW", ID)][, head(.SD, 1), by = "cluster"]
ggData <- ggData[, c("cluster" , "ID")]
ggData$ID <-gsub("HAY_BONE_MARROW_" , "", ggData$ID)
colnames(ggData)[2] <- "HAY_BONE_MARROW"
oupAnnot <- ggData[oupAnnot, on = "cluster"]

## Add in annotation

oupAnnot$annotation <- c("HSC","MPP","MyP1","MyP2","LyP1","LyP2","ProB",
                         "ERP","Ery","MKP","EOBM",
                         "pDC1","pDC2","cDC","Monocyte","InflamMono","Stroma",
                         "PreB","MemoryB","MatureB","BT","PlasmaCell",
                         "CD8+Naive","CD4+Naive","CD4+TCM",
                         "CD8+GZMK1","CD8+GZMK2","CD8+TEMRA","NK")
  
#  rep(c("HSC","MPP","MyP1","MyP2","LyP1","LyP2","ProB"), nrow(oupAnnot))

## output table
png("../Microglia_subtypes_mic_scRNA/findings/08_identifying_marker_genes/clustReAnnotTable.png", width = 12, height = 9, units = "in" 
    , res = 300 )
p1 <-tableGrob(oupAnnot, rows = NULL)
grid.arrange(p1)
dev.off()

### Add new annotation into seurat and replot
tmpMap <- oupAnnot$annotation
names(tmpMap) <- oupAnnot$cluster
newObject$celltype <- tmpMap[newObject@assays$RNA, on = "cluster"]
