library(tidyverse)
library(cowplot)
library(Seurat)
#library(SeuratData)
#remotes::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
#library(Azimuth)
library(patchwork)
library(extrafont)
library(stringr)
library('tidyr')
#devtools::install_github("iaconogi/enrichCM")
library("enrichCellMarkers")
source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")

#Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
#WT_Project <- "WT_data_Microglia"
#AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
#AZT_Project <-  "AZT_data_Microglia"

#WTdata  <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
#AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
#WT_object <-  CreateSeuratObject(counts = WTdata, project = "Microglia subtypes_WT", min.cells = 3, min.features = 300)
#AZT_object <- CreateSeuratObject(counts = AZTdata, project = "Microglia subtypes_AZT", min.cells = 3, min.features = 300)

#integrated_object <-  merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), merge.data= TRUE)

#integrated_object[["strain"]] <- factor(integrated_object@meta.data$orig.ident)
#integrated_object$strain <- str_replace(integrated_object$strain, pattern = "Microglia subtypes_WT", replacement = "WT" )
#integrated_object$strain <- str_replace(integrated_object$strain,  pattern ="Microglia subtypes_AZT", replacement = "AZT" )

#integrated_object[["percent.mt"]] <-PercentageFeatureSet(integrated_object, pattern = "^MT-")

# integrated_object <- subset(integrated_object, subset = nFeature_RNA > 600 & percent.mt < 8)




#GetAssayData(integrated_object@assays$RNA)

#integrated_object[["joind"]] <- Jointlyers(integrated_object, assay = "RNA")

####### ribosomal gene

#ribo.genes <- grep(pattern = "^Rp[s1][[:digit:]]", x = rownames(integrated_object@assays$RNA), value = TRUE)
#integrated_object$percent.ribo <- PercentageFeatureSet(integrated_object,features = ribo.genes)


#some of Microglia genes list ( "lfit3","lfitm3", "lrf7","Hexb" ,"Rps21", "Rps24","Lpl","Apoe","Ccl4","Ccl3", "C3ar1", "Stmn1","Top2a", "Birc5")

#microglia.gene.list <-  c(  "Cx3cr1", "Ctss", "Tmem119", "P2ry12" ,"Cd81" ,"Cst3","Cst7")

#integrated_object$percent.microglia <- PercentageFeatureSet(integrated_object, features = microglia.gene.list)




####### All myeloid cells 
# path= paste("../Microglia_subtypes_mic_scRNA-/findings/01a_QC_strains/", "/",sep = "")

integrated_object <- integrated_object %>% NormalizeData() 
integrated.strain <- subset(integrated_object ) %>% NormalizeData()

 integrated.strain <- JoinLayers(integrated.strain)

 integrated.strain <- AddModuleScore(object = integrated.strain, features = microglia.gene.list)
 
integrated.strain<- integrated.strain %>% NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("batch","ribo.genes", "percent.mt", "nFeature_RNA", "strain","percent.microglia")) %>%
  RunPCA()
integrated.strain <- JackStraw(integrated.strain, num.replicate = 30, dims = 30)
integrated.strain<-ScoreJackStraw(integrated.strain,dims = 1:30)


#######################
ElbowPlot(integrated.strain,ndims = 30) + ggtitle(label = paste("AZT", sep = ""))
ggsave(paste(global_var$global$Path_QC_Strain_findings,"ElbowPlot", ".png", sep = "" ), width = 7, height = 4, dpi = 150)
print(integrated.strain[["pca"]], dims = 1:30, nfeatures = 30)

# Choose dimensions of PCA after checking ElbowPlot and positive genes in each Clusters
pca_dim <- 11
integrated.strain <- integrated.strain %>%
      RunUMAP(reduction = "pca", dims = 1:pca_dim) %>%
      FindNeighbors(reduction = "pca", dims = 1:pca_dim) %>%
      FindClusters(resolution = 0.5)

res="res05" # specify which resolution used to cluster cells

# check dimension reduction

DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 0.001) +
  ggtitle(label = "AZT & WT") + coord_fixed()
ggsave(paste(global_var$global$Path_QC_Strain_findings, "DimPlot1",".png", sep = ""), units = "in", width= 7.3 , height = 7, dpi= 150)


DimPlot(integrated.strain, reduction = "umap" ,label = FALSE, group.by = "strain", pt.size = 0.001 )

DimPlot(integrated.strain, reduction = "umap" ,label = FALSE, group.by = "Strain", pt.size = 0.001 )

ggsave(paste(global_var$global$Path_QC_Strain_findings, res, "_", "DimPlot3", ".png", sep = ""), units = "in", width = 7.3, height = 7, dpi = 150)

# DimPlot(integrated.strain, reduction = "umap", label = TRUE, pt.size = 0.001, split.by = "seurat_clusters")

###### QC: violin plot of nfeatures_RNA, percent.mt, percent.ribo for each cluster
p_QC <-c("nFeature_RNA", "percent.mt", "percent.ribo") %>% map(~QC_plot(integrated.strain@meta.data, .))
p <-plot_grid(plotlist = p_QC, ncol = 1, align = "hv")
title <- ggdraw()+ draw_label(paste(global_var$global$strain,global_var$global$round, global_var$global$res, "QC", sep = " "), fontface = 'bold')
plot_grid(title, p, ncol = 1, rel_heights = c(0.1,1))
ggsave(paste(global_var$global$Path_QC_Strain_findings, global_var$global$strain,"_", global_var$global$round, "_", global_var$global$res, "_", "QC.png", sep = "" ), 
       units = "in" , width = 10, height = 5, dpi = 150 )

##### Find cluster markers 
integrated_joint <- JoinLayers(integrated.strain)

integrated_markers <-FindAllMarkers(integrated_joint, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 300 )
integrated_markers <- integrated_markers %>% rownames_to_column(var = "symbol")

## save cell metadata and marker info into rda
meta <- integrated.strain@meta.data %>% select(-starts_with("ribo_"))

meta_integrated_markers <-merge(meta, y= integrated_markers, add.cell.idec= c("meta", "marker"), Project = "meta_markers"  )

# save(meta, integrated_markers, file = paste(global_var$global$path_data,"Meta_Markers.rda", sep = "") )

##### Do NOT build into function 
### Check cell proportions
prop.table(table(Idents(integrated.strain),integrated.strain$orig.ident), margin= 2)
sum_table <- integrated.strain@meta.data %>% group_by(seurat_clusters) %>%
          summarise(N=n(),
                    med_nCount_RNA=median(nCount_RNA),
                    med_nFeature_RNA=median(nFeature_RNA),
                    med_percent.mt=median(percent.mt),
                    med_percent.ribo=median(percent.ribo),
                    med_percent.migroglia=median(percent.microglia))

### check cell type for some cluster (check selected cluster, don't need to check all)
markers_top <- integrated_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
gene_list <- markers_top %>% filter(cluster == 5) %>% select(symbol) %>% unlist()
results = CMenrich(gene.list = gene_list, species = 'mouse')
DT::datatable(results$enrichments)
results$genes[[1]]


#############
### round 2: just microglia

round = "r3_mg"

# use integrated_strain  result from clustering resolution = 0.6

# use integrated_strain  result from clustering resolusion = 0.6


mg_strain <-subset(integrated.strain, idents = c(0:6, 8:10))
mg.strain <- what_dims(object_type = mg_strain, path = global_var$global$Path_QC_Strain_findings , strain = global_var$global$strain, round = global_var$global$round )
## choose pca dimension
pca_dim <- 22

# try resolution = 0.6
mg.strain <- mg.strain %>%
  RunUMAP(reduction = "pca", dims = 1:pca_dim ) %>%
  FindNeighbors(reduction = "pca", dims = 1:pca_dim) %>%
  FindClusters(resolution = 0.6)



res ="res_06" # specify which resolution used to cluster cells

res ="res_06" # specify which resolution used to cluster cellls

strain <- "AZT_WT"
mg.markers <-  markers(object_type = mg.strain, path = global_var$global$Path_QC_Strain_findings, strain = global_var$global$strain, round = global_var$global$round, 
                     res = res)


# check cell proportions

# check cell proporations
prop.table(table(Idents(mg.strain), mg.strain$strain), margin = 2)  
sum_table <- mg.strain@meta.data %>% group_by(seurat_clusters) %>%
  summarise( N=n(),
             med_nCount_RNA=median(nCount_RNA),
             med_nFeature_RNA=median(nFeature_RNA),
             med_percent.mt=median(percent.mt),
             med_percent.ribo=median(percent.ribo),
             med_percent.microglia=median(percent.microglia))

#### check cell type for some clusters (check selected clusters, don't need to check all)
markers_top<- mg.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
gene_list <- markers_top %>% filter(cluster == 11) %>% select(symbol) %>% unlist()
results= CMenrich(gene.list = gene_list, species = 'mouse')
DT::datatable(results$enrichments)
results$genes[[1]]
