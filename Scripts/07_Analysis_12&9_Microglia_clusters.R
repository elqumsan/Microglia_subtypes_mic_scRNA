
library(ggplot2)
library(gridExtra)
########## Quick analysis for just Microglia cells, that's are in Cluster 12 and 9

integrated_object <- subset(x = integrated_object,idents=c( 7, 9 ))


VlnPlot(integrated_object, feature = microglia.gene.list, pt.size = 0, assay = "RNA", stack = T, flip = T, fill.by = "ident", split.by = "strain",
        group.by = "seurat_clusters")

genes <- c("Cst7",  "Apoe", "Cx3cr1", "Tmem119", "Ifit3", "Ifitm3", "Irf7", "Hexb", "Cd81", "Cst3", "Rplp1", "Rps21", "Rps24" 
           , "C3ar1", "Stmn1", "Top2a", "Birc5" ) 

DotPlot(integrated_object, features = genes) + RotatedAxis() +
  theme(axis.title = element_blank()) +
  scale_y_discrete( labels = function(x) str_wrap(x, width = 20) ) +
  coord_flip() +
  theme(axis.text.y =  element_text( face = "bold.italic"))


DimPlot(integrated_object, group.by = c("orig.ident", "seurat_clusters"))



integrated.meta.stat <- integrated_object@meta.data %>%
  mutate( strain=factor(strain, levels = c("Veh", "AZT")),
          new_clusters = ifelse(seurat_clusters %in% 8:19 , "H", as.character(seurat_clusters)),
          new_clusters=factor(new_clusters, levels = c("0","1","2","3","4","5","6", "7", "H"))) %>%
  group_by(strain, new_clusters) %>%
  summarise(Med_nFeature=median(nFeature_RNA),
            Med_percent_mt= median(percent.mt),
            med_percent.microglia=median(percent.microglia),
            
            N=n()) %>%
  group_by(strain,Percent= N/sum(N)*100)


integratedPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#f0E442", "#0072B2", "#D55E00", "#CC79A7")  


p <- integrated.meta.stat %>%
  ggplot(aes(y=Percent, x= strain, color = strain )) +
  geom_boxplot(outlier.size = 0, alpha= 0.5) +
  geom_point(aes(color = strain), position = position_jitterdodge(), alpha=0.8) +
  scale_color_manual(values = integratedPalette) +
  theme_bw() +
  facet_grid(new_clusters ~ strain, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size= 10),
        strip.text = element_text( face = "bold", size = 12),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom"
  )

##### getting top enriched gene for Cluster 5 and 4

########## Top Markers genes for  Cluster 0
markers_0 <- FindMarkers( integrated.strain, ident.1 = 0, only.pos = TRUE, min.pct= 0.25, test.use = "wilcox", 
                          logfc.threshold = 0.25, ident.use= "seurat_clusters")
markers_0 <- markers_0[markers_0$p_val_adj<0.05,]
markers_0 <- markers_0[order(markers_0$avg_log2FC, decreasing = T),]
(markers_0)[0:10,]
top_gene_names_Cl_0 <-rownames(markers_0[0:10,])
Top_cl_0 <- (top_gene_names_Cl_0)

FeaturePlot(integrated.strain, features = c(Top_cl_0), ncol = 2 )


DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 4, label = FALSE)

DoHeatmap(integrated.strain, features =top_gene_names_Cl_0 , cells = 1:20, 
          size = 6 , angle = 90) + NoLegend()


plot_1 <- FeaturePlot(integrated.strain, features ="Rps5")
HoverLocator(plot = plot_1, information = FetchData(integrated.strain, vars = c("ident", "PC_1", "nFeature_RNA")))


########## Top Markers genes of  Cluster 4
markers_4 <- FindMarkers( integrated.strain, ident.1 = 4, only.pos = TRUE, min.pct= 0.25, test.use = "wilcox", 
                          logfc.threshold = 0.25, ident.use= "seurat_clusters")
markers_4 <- markers_4[markers_4$p_val_adj<0.05,]
markers_4 <- markers_4[order(markers_4$avg_log2FC, decreasing = T),]
head(markers_4)
top_gene_names_Cl_4 <-rownames(markers_4)
Top_cl_4 <- (top_gene_names_Cl_4)[0:10]

FeaturePlot(integrated.strain, features = c(Top_cl_4), ncol = 2 )


DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 4, label = FALSE)

DoHeatmap(integrated.strain, features =Top_cl_4 , cells = 1:20, 
          size = 6 , angle = 90) + NoLegend()


plot_1 <- FeaturePlot(integrated.strain, features ="Rps5")
HoverLocator(plot = plot_1, information = FetchData(integrated.strain, vars = c("ident", "PC_1", "nFeature_RNA")))

####### Top markers genes for cluster 5
markers_5 <- FindMarkers( integrated.strain, ident.1 = 5, only.pos = TRUE, min.pct= 0.25, test.use = "wilcox", logfc.threshold = 0.25, ident.use= "seurat_clusters")
markers_5 <- markers_5[markers_5$p_val_adj<0.05,]
markers_5 <- markers_5[order(markers_5$avg_log2FC, decreasing = T),]
head(markers_5)
top_gene_names_Cl_5 <- rownames(markers_5)
top_10_genes_Cl_5 <- (top_gene_names_Cl_5)[0:20]
FeaturePlot(integrated.strain, features = c(top_10_genes_Cl_5), ncol = 2)

DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 5, label = FALSE)

##VariableFeatures(integrated.strain)[1:100]
DoHeatmap(integrated.strain, features =top_gene_names_Cl_5 , cells = 1:20, 
          size = 6 , angle = 90) + NoLegend()

plot_1 <- FeaturePlot(integrated.strain, features =VariableFeatures(integrated.strain))
HoverLocator(plot = plot_1, information = FetchData(integrated.strain, vars = c("ident", "PC_1", "nFeature_RNA")))

####### Top markers genes for cluster 15
markers_15 <- FindMarkers( integrated.strain, ident.1 = 15, only.pos = TRUE, min.pct= 0.25, test.use = "wilcox", logfc.threshold = 0.25, ident.use= "seurat_clusters")
markers_15 <- markers_15[markers_15$p_val_adj<0.05,]
markers_15 <- markers_15[order(markers_15$avg_log2FC, decreasing = T),]
head(markers_15)
top_gene_names_Cl_15 <- rownames(markers_15)
top_10_genes_Cl_15 <- (top_gene_names_Cl_15)[0:10]
FeaturePlot(integrated.strain, features = c(top_10_genes_Cl_15), ncol = 2)

DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 5, label = FALSE)

##VariableFeatures(integrated.strain)[1:100]
DoHeatmap(integrated.strain, features = top_10_genes_Cl_15, cells = 1:20, 
          size = 6 , angle = 90) + NoLegend()

plot_1 <- FeaturePlot(integrated.strain, features =VariableFeatures(integrated.strain))
HoverLocator(plot = plot_1, information = FetchData(integrated.strain, vars = c("ident", "PC_1", "nFeature_RNA")))



######## Neuron/Neuroblast Tubb3, Meg3, Dcx
FeaturePlot(integrated.strain, features = c("Tubb3", "Meg3", "Dcx"), ncol = 3)

####### Astrocyte/Oligo/Endothelial: Gfap, Olig1, Vtn
FeaturePlot(integrated.strain, features = c("Gfap", "Olig1", "Vtn"), ncol = 3)

####### Macrophage/Monocyte:- F13a1, H2-Aa, Mgl2, Lyve1, and Ccr2
FeaturePlot(integrated.strain, features = c("F13a1", "H2-Aa", "Mgl2", "Lyve1", "Ccr2"), ncol = 3)



##### Percentage of each cluster in Seurat ("making a stacked bar graph of comparison of each condition") 

pt <- table(Idents(integrated.strain), integrated.strain$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)


p <-  ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") + 
#  scale_fill_manual(values = brewer.pal( 21, "Paired")) 
 theme(legend.title = element_blank())

ggsave(p, width = 14, height = 6, filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/precentage_per_cluster.png")


cluster_proportions <- table(integrated_object@active.ident) / ncol(integrated_object)
print(cluster_proportions)
png( filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/cell_proportion_per_cluster.png" )
 barplot(cluster_proportions,
        xlab = "Cluster",
        ylab = "Proportion of cells",
        main = "cell proportion per Cluster"
        )

dev.off()

p <- ggplot(pt, aes(fill= pt$Var2, y= pt$Freq, x= pt$Var1  )) +
  geom_bar(position = "dodge", stat = "identity")
ggsave(p, width = 14 , height = 6 , 
       filename = "../Microglia_subtypes_mic_scRNA/findings/04a_microglia_clustering/cell_proportion_per_cluster_2.png")


 
######### Analysis of transcriptomic signature for each cluster, especially for AZT-related clusters (Gene expression analysis)cluster_markers 
cluster_markers <- FindAllMarkers(integrated.strain, only.pos = TRUE, min.pct = 0.25 , logfc.threshold = 0.25 )
i=0

for (i in levels(cluster_markers$cluster)) {
  cluster_markers[cluster_markers$cluster == i,] <- cluster_markers[cluster_markers$cluster == i,]
 cluster_markers[cluster_markers$p_val_adj,] <-cluster_markers[cluster_markers$p_val_adj< 0.05,]
 cluster_markers[cluster_markers$avg_log2FC,] <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = T),]
  
 #top_gene_names_Cl <- rownames(cluster_markers)[cluster_markers$cluster == i]
 #top_10_genes_Cl <- head(top_gene_names_Cl)
 #FeaturePlot(integrated.strain, features = c(top_10_genes_Cl))
 
 
# DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 5, label = FALSE)
 
 ##VariableFeatures(integrated.strain)[1:100]

 
}

DoHeatmap(integrated.strain, features =top_gene_names_Cl_5 , cells = 1:20, 
           size = 6 , angle = 90) + NoLegend()


####################### Perform differential expression analysis for each cluster

all_markers <- FindAllMarkers(integrated.strain, only.pos = TRUE, min.pct = 0.25 ,logfc.threshold = 0.25)
                         
## Filter for upregulated genes (Positive logFC)
upregulated_genes <- all_markers %>% filter(avg_log2FC > 0 )

upregulated_genes <- all_markers %>% filter(avg_log2FC > 0,  cluster == c(0, 4, 5, 15) )

## Filter for dowregulated genes (negative LogFC)
downregulated_genes  <- all_markers %>% filter(avg_log2FC < 0, cluster == c(0, 4 , 5, 15))

#### Counnt Up/downregulated genes per cluster
cluster_counts <- upregulated_genes %>% group_by(cluster) %>%
  summarise(upregulated_count = n()) %>%
  left_join(downregulated_genes %>% group_by(cluster) %>% summarise(downregulated_count = n()), by = "cluster")

 print(cluster_counts, n = 40)
                         

 ########### Graphically represent for Up/Downregulated genes per clusters   
 ## Plot the data with table Grob
 ss <- tableGrob(cluster_counts)
 
 # Make a scatterplot of your data
 k <- ggplot(cluster_counts) +
   geom_point(  aes(x= cluster_counts$cluster, y=cluster_counts$upregulated_count ))

 ##### Arrange them as you want with grid.arrange 
 grid.arrange(k, ss)

 
 ####### Identification of highly variable genes  
 pbmc <-FindVariableFeatures( integrated.strain, selection.method = "vst", nfeatures = 2000)
 top10_gene <- head(VariableFeatures(pbmc),10 )
plot1 <- VariableFeaturePlot(pbmc) 
plot2 <- LabelPoints(plot = plot1, points = top10_gene, repel = TRUE)
 plot1 +plot2
 