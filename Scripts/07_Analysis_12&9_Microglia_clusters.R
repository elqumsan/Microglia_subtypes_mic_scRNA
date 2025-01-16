
########## Quick analysis for just Microglia cells, that's are in Cluster 12 and 9

integrated_object <- subset(x = integrated_object,idents=c( 9 , 12))


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
markers_5 <- FindMarkers( integrated.strain, ident.1 = 5, only.pos = TRUE, min.pct= 0.25, test.use = "wilcox", logfc.threshold = 0.25, ident.use= "seurat_clusters")
markers_5 <- markers_5[markers_5$p_val_adj<0.05,]
markers_5 <- markers_5[order(markers_5$avg_log2FC, decreasing = T),]
head(markers_5)
top_gene_names_Cl_5 <- rownames(markers_5)
head(top_gene_names_Cl_5)

DimPlot(integrated.strain, reduction = "pca", group.by = "seurat_clusters", order = 5, label = FALSE)



