library(tidyverse)
library(cowplot)
library(Seurat)
library(destiny)
library(SingleCellExperiment)
library(scater)

### load object
DefaultAssay(object.integrated) <- "RNA"
DefaultAssay(integrated.strain) <- "RNA"
DefaultAssay(mg.strain) <- "RNA"

sum_table <-  mg.strain@meta.data %>% group_by(seurat_clusters) %>% 
          summarise( N=n(), ave_nCount_RNA=median(nCount_RNA), ave_nFeature_RNA=median(nFeature_RNA), ave_percent.mt=median(percent.mt), 
                     ave_percent.microglia=median(percent.microglia))
prop.table(table(Idents(mg.strain),mg.strain$strain), margin = 2)


#### Plot both genotypes in all strains (all replicates combined)
# generate meta data, 
integrated.meta <- mg.strain@meta.data %>%
  mutate(strain = factor(strain, levels = c("Veh","AZT")),
         new_clusters= ifelse(seurat_clusters %in% 8:17, "H" ,as.character(seurat_clusters)),
        # new_clusters= ifelse(seurat_clusters %in% 0:0, "H" ,as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("0","1", "2","3","4", "5", "6", "7", "8", "H"))) %>%
  group_by(strain , new_clusters) %>%
  arrange(strain) %>%
  summarise(N=n())

p <- ggplot(integrated.meta, aes(y=N, x=strain , fill= new_clusters )) +
      geom_bar(stat = "identity", position = "fill", color= "black") + 
      labs(y="Fraction", fill = "Clusters") +
      facet_grid(~ strain ) +
    #  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 , face = c("bold", "bold.italic")),
        axis.title.x  = element_blank(),
        strip.text.x =  element_text(face = "bold"),
        axis.ticks.x = element_blank(),
        axis.line.x =  element_blank()
  #      legend.position = "bottom"
      )
ggsave(paste(global_var$global$path_microglai_statistics, "fraction_replicates_seperated.png", sep = "/"), p , width = 3.5 , height = 5 , units = "in")

############ Box Plot for all microglia
############ generate meta data, for statistical testing and box plot 
integrated.meta.stat <- mg.strain@meta.data %>%
                mutate( strain=factor(strain, levels = c("Veh", "AZT")),
                        new_clusters = ifelse(seurat_clusters %in% 8:17 , "H", as.character(seurat_clusters)),
                        new_clusters=factor(new_clusters, levels = c("0","1","2","3","4","5","6", "7", "H"))) %>%
                group_by(strain, new_clusters) %>%
                summarise(Med_nFeature=median(nFeature_RNA),
                          Med_percent_mt= median(percent.mt),
                          med_percent.microglia=median(percent.microglia),
                          
                          N=n()) %>%
  group_by(Percent= N/sum(N)*100)

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

ggsave(paste(global_var$global$path_microglai_statistics, "cluster_box_all.png", sep = "/"), p , width = 4, height = 10, units = "in")


####### Perform two-way ANOVA to determine the effect of strain on the percent of microglia subclusters

clusters <- unique(integrated.meta.stat$new_clusters) %>% as.list()
data  = integrated.meta.stat %>%  filter(new_clusters %in% clusters[1])
aov_object <-aov( Percent ~ strain, data = data)
aov.pvals <- summary(aov_object)


aov.pvals = aov.pvals[[1]][[2]] %>% t() %>% as.data.frame()
names(aov.pvals)<- c( "Residuals")
aov.pvals <- aov.pvals %>% 
              select(-Residuals) %>%
           mutate(Cluster = clusters[1] %>% as.character())

aov_Strain <- function(cluster, data){
         data = data %>% filter(new_clusters %in% clusters) 
        aov_object = aov(Percent ~ strain, data= data)
        aov.pvals = summary(aov_object)
        aov.pvals = aov.pval[[1]][[2]] %>% t() %>% as.data.frame()
        names(aov.pvals) <- c("Residuals")
        aov.pvals <- aov.pvals %>% 
          select(-Residuals) %>%
          mutate(Cluster = cluster %>% as.character())
        return(aov.pvals)
}


aov_Strain_object <- function(cluster, data){
  data = data %>% filter(new_clusters %in% cluster)
  aov_object = aov(Percent ~ strain , data= data)
  return(aov_object)
  
}

aov_strain_table <- clusters %>% map_df(aov_Strain, data = integrated.meta.stat )

######### check the statistics on nFeature, percent of microglia, and percent of ribosomal genes for each cluster 

#### nFeature
aov_stat <- aov(Med_nFeature ~ new_clusters, data = integrated.meta.stat)
aov_table <-TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>%
  mutate(comparison = paste(" ", comparison , sep = ""), Significance=ifelse(p.adj < 0.05, "S", "NS"))
write_delim(aov_table, path = paste(global_var$global$path_microglai_statistics, "Med_nFeature_comp_cluster.txt", sep = "/"), delim = "\t")
filter(aov_table, Significance == "S")


###### percent of mitochodria
aov_stat = aov(Med_percent_mt ~ new_clusters, data = integrated.meta.stat)
aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>%
  mutate(comparison = paste(" ", comparison, sep = " "), Significance=ifelse(p.adj < 0.05 , "S", "NS"))
write_delim(aov_table, path = paste(global_var$global$path_microglai_statistics, "Med__percent_mitochondria.txt", sep = "/"), delim = "\t")
filter(aov_table, Significance =="S")

#write_delim(aov_table, path = paste(global_var$global$path_microglai_statistics, "Med_percent_mt_comp_cluster.txt", sep = "/"), delim = "\t")

################## percent of Microglia
aov_stat = aov(med_percent.microglia ~ new_clusters, data = integrated.meta.stat)

aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var =  "comparison") %>%
  mutate(comparison = paste(" ", comparison, sep = " "), Significance=ifelse(p.adj < 0.05 , "S", "NS"))
write_delim(aov_table, path = paste( global_var$global$path_microglai_statistics, "med_percent.microglia.txt", sep = "/"), delim= "\t")

#######################   Pseudotime analysis (diffusion map ) 
##################### too many cells for diffusion map, need sampling 

integrated.strain$final_clusters <-ifelse(integrated.strain$seurat_clusters %in% 0:0, "H", 
       integrated.strain$seurat_clusters %>% as.character())

sampling <- integrated.strain@meta.data %>% 
              rownames_to_column(var = "cell_ID") %>%
              group_by(strain) %>%
              sample_n(1000)  # take 1000 random cells from each group

mg.small <- subset(integrated.strain, cells=sampling$cell_ID)

mg.small <- as.SingleCellExperiment(mg.small)

# Use diffusion map to calculate pseudotime
pca <- reducedDim(mg.small)
cellLables <- mg.small$seurat_clusters

pca_tidy <- as.data.frame(pca) %>% rownames_to_column()

rownames(pca)<- cellLables

dm <- DiffusionMap(pca)

dpt <- DPT(dm)

mg.small$pseudotime_dpt <- rank(dpt$dpt)

df <- colData(mg.small) %>% as.data.frame()

df$final_clusters <- ifelse(df$seurat_clusters %in% 0:0, "H", df$seurat_clusters %>% as.character())


ggplot(df, aes(pseudotime_dpt, fill= final_clusters)) +
  geom_histogram(binwidth = 100 , color = "grey", size=0.1) +
  facet_grid(strain ~., switch = "y") +
  scale_y_continuous("count", position = "right") +
  labs(x="DAM <- pseudotime -> Homostatic ") +
  theme_bw() +
  theme(text = element_text(family = "Arial", size = 10),
        strip.text.y = element_text(size = 5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y.right = element_blank()
       
        )
ggsave(paste(global_var$global$path_microglai_statistics, "pseudotime.png", sep = "/"), width = 3.5 , height = 5.3 , units = "in", dpi = 600)
