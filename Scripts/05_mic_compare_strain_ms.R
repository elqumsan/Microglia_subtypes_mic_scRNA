####
####
#### Comparing cluster marker genes across strains with reported microglia 
library(tidyverse)
library(UpSetR)
library(Seurat)
library(cowplot)
library(devtools)
library(Vennerable)
library(extrafont)



#### set up marker gene sets for each cluster 

cluster_name <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")

strain = c("Veh", "AZT")


##### Microglia markers 
genes <- c("Tmem119", "Cx3cr1", "P2ry12", "Sall1", "Gpr34", "Spi1", "Mafb" , "Maf", "Mef2a", "Irf8", "Fcrls" )
Idents(integrated.strain) <- "Strain"
p <- genes %>% 
  map(~VlnPlot(integrated_object %>% subset(subset = seurat_clusters %in% c( "8" , "9", "10", "11", "12")), features =. , pt.size = 0 , split.by= "strain" , ncol = 4 
                   , cols = c("#E69F00", "#999999") ) +
                     theme( legend.position = "none",
                            axis.title = element_blank(),
                            axis.text.x =  element_text(size = 8),
                            axis.text.y = element_text(size = 10),
                            title = element_text( size = 10, family = "Arial")
                     ) )

                   
plot_grid(p[[1]], p[[2]], p[[3]],p[[4]], p[[5]],p[[6]], p[[7]],p[[8]], p[[9]], p[[10]], p[[11]], 
          align = ("hv"), nrow = 3, ncol = 4)

ggsave(paste(global_var$global$path_compare_strain,"/", "microglia_markers.png", sep = "" ), 
       units = "in" , width = 6.4, height = 5, dpi = 300 )


### Z-score plotting using split violin and grouped by cluster 
inpalette  <- c("#999999", "#E69F00", "#56B4E9", "#009E73" ,"#F0E442", "#0072B2", "#D55E00", "#CC79A7")
comp <- "homeo1"
cluster <- c("6" ,"7", "8", "9", "10", "11", "12","13" ,"14","15")
strain <- c("Veh", "AZT") %>% rev()
width_height <- c(3,5)

integrated.strain@meta.data %>% 
   filter(seurat_clusters %in% cluster  ) %>%
  
  mutate( suerat_clusters = factor(seurat_clusters, levels =  c("6" ,"7","8", "9","10" ) ),
    strain = factor(strain, levels = c("Veh","AZT")) ) %>%

 ggplot(aes(y= "", x=strain, fill= strain)) + #change y to the same comp but no quotation mark
     facet_grid( seurat_clusters ~ . , scales = "fixed") +
 introdataviz::geom_split_violin( trim = TRUE) +
  geom_boxplot(width = 0.25, notch =  FALSE, notchwidth = .4 , outlier.shape = NA , coef= 0) +
  scale_fill_manual(values = inpalette) +
  
  coord_flip() +
  theme_bw() +
  labs(x=NULL,  y = "z-score") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12, family = "Arial"),
        axis.text = element_text(size = 10),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold", size = 12)
        )

ggsave(paste(global_var$global$path_compare_strain, "enrichment_Z_score.png", sep="/"), dpi = 300, width = 6.4, height = 5)

#############
## microglia markers

genes <- c("Tmem119", "Cx3cr1", "P2ry12", "Sall1", "Gpr34", "Spi1",  "Mafb", "Maf", "Mef2a", "Irf8",  "Fcrls", "Olfml3")
Idents(integrated.strain) <- "Strain"

p<- genes %>% 
  map(~VlnPlot(integrated.strain %>% subset(subset = seurat_clusters %in% c("3", "4","5","6", "7", "8")) , features = . , pt.size = 0 , 
               split.by = "strain" , ncol = 4, cols = c("#E69F00", "#999999")) +
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_text(size = 10),
              title = element_text(size = 10, family = "Arial")
              )
        )

plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]] , p[[10]], p[[11]], p[[12]],
          align = ("hv"), nrow = 3, ncol = 4)


######################################

meta <- integrated.strain@meta.data %>% select(-starts_with("ribo_"))
head(integrated_markers)
sample_IDs<- colnames(integrated.strain)
 
final_meta  <- data_frame(meta[1:25167,],integrated_markers[1:25167,], sample_IDs)

head(final_meta)



#### set up marker gene sets for each cluster
cluster_name <- c("Cluster1", "Cluster20", "Cluster21", "Cluster22")
strain= c("WT", "AZT")
marker_genes<-  FindAllMarkers(integrated_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 300)
marker_genes <- marker_genes %>% rownames_to_column(var = "symbol")
marker_genes <- marker_genes %>% group_by(cluster) %>% top_n(n=10, wt= p_val_adj) # do not put quotation mark on the wt argument #  sho top markers 

### Wrap into a function to retrieve DE gene sets

### Wrap into a function to retrieve DE gene set

cluster_DE <- function(logFDR_cut, logFC_cut, df){
  # input:
  # LogFDR_cut, threshold of -log10FDR
  # logFC_cut, threshold of log2FC
  # df: DE gene table for each cluster "marker_genes"
  # output:
  # a list containing marker genes for a given cluster for all strains
  # strain: global environment
  DE_list <- vector(mode = "list", length = length(summed))

  # df: DE gene table for each cluster 
  # output:
  # a list containing marker genes for a given cluster for all strains
  # strain: global environment
  DE_list <- vector(mode = "list", length = length(strain))

  
  if(logFC_cut > 0){
    for(i in seq_along(strain)){
      DE_list[[i]] <- marker_genes %>% select(symbol, contains(strain[i]))  %>%

      DE_list[[i]] <- df %>% select(symbol, contains(strain[i]))  %>%
        filter( str_detect(symbol, "^Gm", negate = TRUE) ) %>%
        filter_at(var(contains("p_val_adj")), any_vars(-log10(.)> logFDR_cut )) %>%
        filter_at(vars(contains("logFC")), any_vars(. > logFC_cut)) %>%
        select(symbol) %>% units()
    }
  } else {
    for ( i  in seq_along(strain)) {
      DE_list[[i]] <- marker_genes %>% select(symbol , contains(strain[i])) %>%

      DE_list[[i]] <- df %>% select(symbol , contains(strain[i])) %>%
       filter(str_detect(symbol, "^Gm", negate = TRUE)) %>%
        filter_at(vars(contains("p_val_adj")), any_vars(-log10(.)> logFDR_cut)) %>%
        filter_at(vars(contains("logFC")), any_vars(. < longFC_cut)) %>%
        select(symbol) %>% unlist()
    }
        
      }
  names(DE_list) <- strain
  return(DE_list)

}  # End of cluster_DE function 


cluster_DE(17 , 0.25, 2)

cluster_DE(18, 0.25,2)
