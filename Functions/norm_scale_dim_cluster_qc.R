



what_dims <- function(object_type , path, strain, round){
       
  object_type <- JoinLayers(object_type)
  object_type <-AddModuleScore(object = object_type , features = microglia.gene.list, ctrl = 100, name = 'Microglia_Features' )
  object_type <- object_type %>%
    NormalizeData()%>%
    FindVariableFeatures(selection.methods = "vst", nFeatures = 3000) %>%
    ScaleData(vars.to.regress = c("percent.ribo", "percent.microglia","percent.mt", "nFeature_RNA")) %>%
    RunPCA()
  
  object_type <- JackStraw(object_type, num.replicate = 30, dim = 30)
  object_type <-ScoreJackStraw(object_type, dims = 1:30)
  ElbowPlot(object_type, ndims = 30) + ggtitle( label = paste(strain, round , sep = ""))
  ggsave(paste(path, strain,"_", round, "_", "ElbowPlot.png", sep = "" ), units = "in", width = 7, height = 4, dpi = 150)
  
  print(object_type[["pca"]], dims = 1:30 , nfeature = 30)
  
  return(object_type)
}


# QC: Violin plot of nfeatures_RNA, percentmt, for ecah cluster

QC_plot <- function(data, y ){
  
  p <- ggplot(data, aes_string(x= "seurat_clusters", y= y  , colors= "seurat_clusters")) + 
    geom_violin()+
    geom_boxplot(width= 0.07, outlier.shape = NA, color = "red", alpha= 0.7) +
    theme_bw()+
    theme(legend.position = "none", axis.title.x =  element_blank())
  return(p)
  
}

##### Find markers 

markers <- function(object_type, path, strain, round, res ){
  # check dimension reduction 
  DimPlot(object_type, reduction = "umap", label = TRUE, pt.size = 0.001) +
    ggtitle(label = strain ) + coord_fixed()
  ggsave(paste(path, strain, "_", round , res , "_", "DimPlot1.png", sep = "" ), units = "in", width = 7.3, height = 7, dpi = 150 )
  
  ##### QC:violin plot of nfeature_RNA, percent.mt, percent.ribo for each cluster
  ## function: QC_plot(data, y) from source
  p_QC <- c("nFeature_RNA", "percent.mt", "percent.ribo") %>% map(~QC_plot(object_type@meta.data, .))
  p <- plot_grid(plotlist =p_QC,ncol = 1, align = "hv")
  title <- ggdraw() + draw_label(paste(strain, round, res, "QC", sep = ""), fontface = 'bold')
  plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
  ggsave(paste(path, strain, "_", round, "_", res, "_", "QC.png", sep = "" ), units = "in", width = 10, height = 5 , dpi = 150)
  
  # Find cluster markers
  object_markers <- FindAllMarkers(object_type, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25 , max.cells.per.ident = 300) # max cells.per.ident
  object_markers <- object_markers %>% rownames_to_column(var = "symbol")
  
  # save cell metadata and marker info into rda
  meta <- object_type@meta.data %>% select(-starts_with("ribo_"))
  save(meta, object_markers, file = paste(path, strain, "_", round, "_", res , "_", "Meta_Marker.rda", sep = ""))
  return(object_markers)
  
}
  