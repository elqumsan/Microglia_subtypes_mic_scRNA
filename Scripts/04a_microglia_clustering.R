
library(tidyverse)
library(cowplot)
library(Seurat)
library(patchwork)


integrated_object@meta.data %>% group_by(strain) %>% summarise(N=n())
sum_table <- integrated_object@meta.data %>% group_by(seurat_clusters) %>% summarise(N= n(), ave_nCount_RNA=median(nCount_RNA) ,
                           ave_nFeature_RNA=median(nFeature_RNA), ave_percent.mt=median(percent.mt), ave_percent.microglia=median(percent.microglia)  )
prop.table(table(Idents(integrated_object), integrated_object$strain), margin = 2)

Idents(integrated.strain) <- "seurat_clusters"
plot_title = "mic_cluster"
i=17 #(PCA dim)
j= 0.6
j=0.6 #(resolution)

##### UMAP Plot
### UMAP Plot - Strain integrated

DimPlot(integrated_object, reduction = "umap", label = TRUE, pt.size = 0.001, label.size = 5) +
  ggtitle( label = paste("integration", j, "res", sep = "_"))
  coord_fixed() +
  theme(axis.title = element_blank(), legend.position = "none")
ggsave(paste(global_var$global$path_microglia_clustering, "/" , "umap_", i, "_res_", j, "_Dimplot_Strain_small", ".png", sep=""), units = "in", 
       width = 4 , height = 4, dpi = 300 )


##### check marker genes of each cluster, list top 10 marker genes 

integrated_markers %>% group_by(cluster) %>% top_n(-10, wt = p_val_adj)


####### UMAP plot strain  spitted 

# Idents(integrated.strain) <- "final_clusters"

p <- DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 1E-10, split.by = "strain", ncol = 2, label.size = 3 , repel = TRUE ) +
  ggtitle(label = paste("AZT agnaist WT Clustering")) +
  coord_fixed() +
  theme(title = element_text(size = 18),
        axis.text =  element_text(size = 12),
        legend.position = "none")

ggsave(paste(global_var$global$path_microglia_clustering, "/", "Dimplot_integarted_data.png", sep = "" ), p , dpi = 300)



######### Dot Plot of microglia marker genes for each clusters 

#Idents(integrated.strain) <- "final_clusters"
 genes <- c("Fcrls","P2ry12","Cx3cr1", "Trem2", "C1qa" ,"Tmem119" ) 
 
# genes <- c("Cx3cr1", "Tmem119", "Hexb", "Cd81", "Cst3", "Rplp1", "Rps21", "Rps24", "Apoe" , "Nos2", "Tspo",  "Cd86" , "Cd16", "Cd32", 
#            "Cd40", "IL6" , "IL12", "IL18", "Mmp9", "Mmp2", "Tnfrsf1a", "Ccl5", "Ccl10", "Cxcl1", "Cxcl9","Cxcl10",
#            "Arg1", "IL4", "IL10", "IL13", "Ccl2", "Ccl17", "Ccl22", "Ccl24", "Ym1", "Fizz1", " Cd163", "Cd206", "Tgfb1",)
 
 file_name <- paste(global_var$global$path_microglia_clustering, "/Dotplot_all_gene.png", sep = "")
 DotPlot(integrated.strain, features = genes) + RotatedAxis() +
   theme(axis.title = element_blank()) +
   scale_y_discrete( labels = function(x) str_wrap(x, width = 20) ) +
   coord_flip() +
   theme(axis.text.y =  element_text( face = "bold.italic"))

 ggsave( file_name,  dpi = 300) 
 
 ####### Feature plots of marker genes of microglia subclusters 
 
 genes <- c("Fcrls","P2ry12","Cx3cr1", "Trem2", "C1qa" ,"Tmem119")

 ### Plot one gene
p  <- genes %>% 
   map(~FeaturePlot(integrated.strain, features = ., min.cutoff = "q9" ,label = TRUE, repel = FALSE, ncol = 2, order = TRUE)+
         coord_fixed() +
         theme(axis.line = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank())  +
         scale_y_discrete(labels= function(x) str_wrap(x, width = 20)) 
         
       )
 
(p[[1]]+p[[2]])/(p[[3]] +p[[4]])/(p[[5]] + p[[6]])

ggsave(paste(global_var$global$path_microglia_clustering, "/Feature_plot_all_2.png", sep = ""), units = "in", width = 10, height = 7 , dpi = 300)


### Check immediate early genes (IEG)
### To verify that our microglia prepared by mechniacal dissociation 

### Feature plots for IEG
genes <- c("Fos", "Fosb", "Dusp1", "Nr4a1" , "Arc", "Egr1")

p <- genes %>% map(~FeaturePlot(integrated.strain, features =., min.cutoff = "q9", label=TRUE, repel=TRUE, ncol= 2, order= FALSE)+
                coord_fixed()+
                theme(axis.line = element_blank(),
                      axis.title = element_blank(),
                      axis.ticks = element_blank() 
                      ) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 20) )
                  
                )

(p[[1]]+p[[2]])/(p[[3]]+p[[4]])/(p[[5]]+p[[6]])

ggsave(paste(global_var$global$path_microglia_clustering, "/Feature_plot_all_unordered.png", sep = ""), units = "in", width = 10 , height = 7, dpi = 300)


#### Dot plot for IEG
DotPlot(integrated.strain, features = genes) + RotatedAxis() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text( face = "bold.italic") )+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  coord_flip()
ggsave(paste(global_var$global$path_microglia_clustering, "/Dotplot_IEG.png", sep = ""), units = "in", width = 6.1, height = 3.4, dpi = 300 )


###### Check the gene number, percent of percent of mitochodrial and microglia genes in each cluster 


QC_plot_single2 <- function(data , y ){
  p <- ggplot(data, aes_string( x= "seurat_clusters" , y=y, color= "seurat_clusters"))+
    geom_violin()+
    geom_boxplot(width =0.07, outlier.shape = NA, color = "black", alpha=0.7) +
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.ticks = element_blank())
    return(p)
}

p_QC <- c("nFeature_RNA", "percent.mt", "percent.microglia") %>% map(~QC_plot_single2(integrated.strain@meta.data, .))

p <- plot_grid(plotlist = p_QC, nrow = 3, ncol = 1 , align = "hv")
plot_grid(p, nrow = 1, rel_heights = c(0.1, 1))

ggsave(paste(global_var$global$path_microglia_clustering, "/all_microglia_integrated2.png", sep = ""), units = "in", width = 4.5, height = 3.5, dpi = 300)


VlnPlot(integrated_object, feature = microglia.gene.list, pt.size = 0, assay = "RNA", stack = T, flip = T, fill.by = "ident", split.by = "strain",
        group.by = "seurat_clusters")

ggsave(paste(global_var$global$path_microglia_clustering, "/percent_of_microG_in_each_cluster.png", sep = ""), units = "in", width = 4.5, height = 3.5 ,
       dpi = 300)

