#### Integration Microglia

library(tidyverse)
library(cowplot)
library(Seurat)

#source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")
# meta <- load("../Microglia_subtypes_mic_scRNA-/data/Meta_Markers.rda")


### collect cells from different strains analysis results
rda_file <- c("WT", "AZT")
rda_list <- vector(mode = "list", length = length(rda_file))


## save cell metadata and marker information into rda

## save cell metadata and marker infor into rda

meta <- object.significance.PCA.scores@meta.data %>% select(-starts_with("^Mico_"))


for(i in seq_along(rda_list)) {

  (rda_file[[i]])
  rda_list[[i]] = meta %>% rownames_to_column(var = "cell_id")

  rda_file[[i]]
  rda_list[[i]] <- meta %>% rownames_to_column(var = "cell_id")

  
}

names(rda_list) <- c("WT", "AZT")
rda_list %>% map(dim)


## collect cell ID:-
cells <- rda_list %>% map_df(~ select(., cell_id )) %>% unlist() %>% unname()

# check the number of the cells

rda_list %>% map(dim) %>% map(~ .[1]) %>% unlist() %>% sum() # 50334

## Subset cells from integrated_strain, split by strain, then normalize and find variable genes before integration 

integrated.merged <- integrated.strain %>% 
  subset(cells = cells ) %>%
  SplitObject(split.by = "strain")


## Now integrated.merged is splitted object:

## Now integrated_strain is splitted object:
integrated.merged %>% map(dim)


for (i in seq_along(integrated.merged)) {
  integrated.merged[[i]] <- integrated.merged[[i]] %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
  
}


# integration


integrated.anchors <- FindIntegrationAnchors(object.list = integrated.merged,
                                              dims = 1:30,
                                              anchor.features = 3000)

object.integrated <-IntegrateData(anchorset = integrated.anchors   , dims = 1:30 )


DefaultAssay(object.integrated) <- "integrated"

## Run the standard workflow for visualization and clustering 

object.integrated <- ScaleData(object.integrated, vars.to.regress = c("strain", "percent.mt", "percent.microglia"))
object.integrated <- RunPCA(object.integrated, npcs = 30)

## Cluster cells by first using PCA dimension of 30, later tune the different parameters to find best clustering 
object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:30)
object.integrated <- FindNeighbors(object.integrated, reduction = "pca", dims = 1:30 )
#object.integrated <- FindClusters(object.integrated, resolution = 0.5)

object.integrated <- FindClusters(object.integrated, resolution = c(0.5 ,0.6, 0.7, 0.8 ))


#### determine cell clusters by exploring different PCA dimensions included and cluster resolution


## Use different PCA and resolution to determine clusters less PCA dimension : pca = 25 (previously =30 )
# generate ElbowPlot to have a guess of how many pca take

object.integrated <- JackStraw(object.integrated, num.replicate = 15, dims = 30)
object.integrated <- ScoreJackStraw(object.integrated)

ElbowPlot(object.integrated, ndims = 30)
ElbowPlot(object.integrated, ndims = 30) + ggtitle(label = "Integrated_data")

ggsave(paste(global_var$global$path_microglia_integration, "ElbowPlot_integrated_data", ".png", sep = ""), width = 7, height = 4, dpi = 200)

print(object.integrated[["pca"]], dims = 1:30 , nfeatures = 10)

#integrated.anchors <- FindIntegrationAnchors(object.list = integrated.merged,
#                                              dims = 1:30,
#                                              anchor.features = 3000)

#integrated.merged <-IntegrateData(anchorset =  , dims = 1:30 )





# Use different PCA and resolution to determine clusters less PCA dimension: pca dim 12

#DefaultAssay(object.integrated) <- "integrated"

# generate ElbowPlot to have a guess of how many pca to take
#object.integrated <- JackStraw(object.integrated, num.replicate = 15, dims = 30)
#object.integrated <-ScoreJackStraw(object.integrated, dims = 1:30)
#ElbowPlot(object.integrated, ndims = 30) + ggtitle(label = "Integrated_data")

# print(object.integrated[["pca"]], dims = 1:30 , nfeatures = 10)

# then 
pca_dim <- c(17)
res <- c(0.5, 0.6, 0.7, 0.8)

## Use QC_Plot function in the for Loop below

for(i in pca_dim){
  DefaultAssay(object.integrated) <- "integrated"
  object.integrated <- object.integrated %>%
  RunUMAP(reduction = "pca", dims = 1:i) 
  
 # for(j in res){
    
#    object.integrated <- object.integrated %>%

#  DefaultAssay(object.integrated) <- "integrated"
#  object.integrated <- object.integrated %>%
#  RunUMAP(reduction = "pca", dims = 1:i) 
  
  for(j in res){
    
    object.integrated <- object.integrated %>%
      FindNeighbors(reduction = "pca", dims = 1:i) %>%
      FindClusters(resolution = j)
    
    ## markers  ( only calculate under the lowest resolution )
    
    if(j== 0.5){
      DefaultAssay(object.integrated) <- "RNA"
      ##### Find cluster markers 
      object.integrated<- JoinLayers(object.integrated) ### Have always to be before FindAllMarkers Function 
      
      integrated_markers <- FindAllMarkers(object.integrated, only.pos = TRUE , min.pct= 0.25, logfc.threshold = 0.25 , max.cells.per.ident = 300 ) # max.cells.per.ident
      proptb1 <- prop.table(table(Idents(object.integrated), object.integrated$strain), margin = 2)
 #     save(integrated_markers, proptb1, file = paste(global_var$global$path_microglia_integration,"pcs_", i, "_res_", j, ".rda", sep = ""))
      
    } else {
      proptb1 <- prop.table(table(Idents(object.integrated), object.integrated$strain), margin = 2)
 #     save(proptb1, file = paste(global_var$global$path_microglia_integration, "pca_", i, "_res_", j, "_proptb1", ".rda", sep = ""))
      
    }  # end of else
    
    
    # visualize
    DimPlot(object.integrated, reduction = "umap", label = TRUE, pt.size = 0.001) +
      ggtitle(label = paste("integration" , i, "res" , j, sep = "_")) + coord_fixed()
    ggsave(paste(global_var$global$path_microglia_integration, "umap_",  i, "_res_", j,  ".png", sep = "" ), units = "in", width = 7.3, height = 7 , dpi = 150)
    
    DimPlot(object.integrated, reduction = "pca", label = TRUE, pt.size = 0.001) +
      ggtitle(label = paste("integration" , i, "res" , j, sep = "_")) + coord_fixed()
    ggsave(paste(global_var$global$path_microglia_integration, "pca_",  i, "_res_", j,  ".png", sep = "" ), units = "in", width = 7.3, height = 7 , dpi = 150)
    
    
    p_QC <- c("nFeature_RNA", "percent.mt", "percent.ribo", "percent.microglia"  ) %>% map(~QC_plot(object.integrated@meta.data, .))
    p <- plot_grid(plotlist = p_QC, ncol = 1, align = "hv")
    title <- ggdraw() + draw_label("integration_data" , fontface = 'bold')
    plot_grid(title, p , ncol = 1, rel_heights = c(0.1, 1 ) )
    ggsave(paste(global_var$global$path_microglia_integration, "pca_", i, "_res_", j, "_QC" ,".png", sep = "" ), units = "in", width = 10, height = 5, dpi = 200)
    
  }# end for res loop
  
} # end for pca_dim loop 





#### Then determine the best pca dim and res combination:  run it again to get marker genes

## !! Best pca dimension : 20, 21 , 22

# top_30 <- integrated_markers %>% group_by(cluster) %>% top_n(n= 30, wt= avg_log2FC)

## run the best possible combination
#i=21
#j=0.5

#DefaultAssay(object.integrated) <- "integrated"
#object.integrated <- object.integrated %>%
#  RunUMAP(reduction = "pca", dims = 1:i) %>%
#  FindNeighbors(reduction = "pca", dims = 1:i) %>%
#  FindClusters(resolution = j )

#### markers (only calculate under the lowest resolution )
#DefaultAssay(object.integrated) <- "RNA"

#object.integrated<- JoinLayers(object.integrated)

#integrated_markers <- FindAllMarkers(object.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 300) ## max.cell.per.ident
#proptb1 <- prop.table(table(Idents(object.integrated), object.integrated$strain), margin = 2)

#save(integrated_markers, proptb1, file=paste(global_var$global$path_data, "pca_", i, "_res_" , j,"_markers_proptb1", ".rda", sep = ""))


## save the object with PCA
##saveRDS(object.integrated, paste( global_var$global$path_data, "object.integrated.rds", sep = ""))

#meta <- object.integrated@meta.data %>% rownames_to_column(var="cells")

#saveRDS(meta, paste(global_var$global$path_data, "mate.rds", sep = ""))

#      integrated_markers <- FindAllMarkers(object.integrated, only.pos = TRUE , min.pct= 0.25, logfc.threshold = 0.25 , max.cells.per.ident = 300 ) # max.cells.per.ident
#      proptb1 <- prop.table(table(Idents(object.integrated), object.integrated$))
      
#    } # end of if J
    
#  }# end for res loop
  
#} # end for pca_dim loop  
