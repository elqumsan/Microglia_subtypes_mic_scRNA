
library(tidyverse)
library(recount)
library(SummarizedExperiment)
library(S4Vectors)

library(stringr)
library('tidyr')


library(shiny)
library(bslib)
library(ggExtra)
library(purrr)

library(knitr)
library(SummarizedExperiment)
library(SeuratObject)
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork) 
library("BiocManager")
library("glmGamPoi")
library('enrichCellMarkers')

source("../Microglia_subtypes_mic_scRNA/Functions/norm_scale_dim_cluster_qc.R")

#source("../Microglia_subtypes_mic_scRNA-/Functions/norm_scale_dim_cluster_qc.R")
#source("../Microglia_subtypes_mic_scRNA-/scripts/02_QC_starin_split.R")



Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
WT_Project <- "WT_data_Microglia"
AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
AZT_Project <-  "AZT_data_Microglia"

WTdata  <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
WT_object <-  CreateSeuratObject(counts = WTdata, project = "Microglia subtypes_WT", min.cells = 3, min.features = 300)
AZT_object <- CreateSeuratObject(counts = AZTdata, project = "Microglia subtypes_AZT", min.cells = 3, min.features = 300)

integrated_object <-  merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), merge.data= TRUE)
# object.ref <- subset(integrated_object, orig.ident %in% c("Microglia subtypes_WT" , "Microglia subtypes_AZT"))

strain <- c("Microglia subtypes_WT", "Microglia subtypes_AZT")
cols =  c("#888888", "#00AA00")

integrated_object[["strain"]] <- factor(integrated_object@meta.data$orig.ident, levels = strain)
integrated_object$strain <- str_replace(integrated_object$strain,pattern = "Microglia subtypes_WT", replacement = "WT" )
integrated_object$strain <- str_replace(integrated_object$strain, pattern ="Microglia subtypes_AZT", replacement = "AZT" )

#sample_ID <- colnames(integrated_object)
#meta_data <- data_frame(Sample_ID = sample_ID,Cust_ID ,Strain = integrated_object$Strain )


meta <- integrated_object@meta.data

meta_tidy <- meta %>% 
  select(orig.ident,nCount_RNA, nFeature_RNA, strain) %>%
  gather(-orig.ident, -strain ,key= "QC", value="value" )
#VlnPlot(object.ref, features = c("nFeature_RNA","nCount_RNA" ) )

#map(~QC_plot(integrated_object@meta.data, .))

######### ribosomal gene
ribo.genes <- grep( pattern = "^Rp[s1][[:digit:]]", x = rownames(integrated_object@assays$RNA),value = TRUE)
integrated_object$percent.ribo <- PercentageFeatureSet(integrated_object, features = ribo.genes)
microglia.gene.list <-  c(  "Cx3cr1", "Ctss", "Tmem119", "P2ry12" ,"Cd81" ,"Cst3","Cst7") 
integrated_object$percent.microglia <- PercentageFeatureSet(integrated_object, features = microglia.gene.list )
integrated_object[["percent.mt"]] <- PercentageFeatureSet(integrated_object, pattern = "^MT-")


## Check cell type for some cluster (check selected cluster, don't need to check all)

results = CMenrich(gene.list = c('Sall1', 'Hexb', 'Fcrls', 'Gpr43', 'Cx3cr1', 'Tmem119', 'Trem2', 'P2ry12', 'Mertk', 'Pros1','Siglech'), species = 'mouse' )


##################################
pca_dim = 11
integrated_object<- integrated_object %>% 
                      NormalizeData() %>%
                      FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
                      ScaleData(vars.to.regress = c("nFeature_RNA", "strain","percent.microglia" ,"percent.ribo", "percent.mt")) %>%
                      RunPCA() %>%
                      RunUMAP(reduction = "pca", dims = 1:pca_dim ) %>%
                      FindNeighbors(reduction = "pca", dims = 1:pca_dim ) %>%
                      FindClusters(resolution = 0.5)


#### integrated.strain object has no layers, all layers have been merged to be eligible for the rest analysis 
integrated.strain <- integrated_object
integrated.strain <- JoinLayers(integrated.strain)



#######

mg.strain <-what_dims(object_type = integrated_object, path = global_var$global$Path_QC_Strain_findings, strain = global_var$global$strain , round = global_var$global$round  )

cells <- WhichCells(integrated_object)

CellsMeta = integrated_object@meta.data
randomnumbers <- runif(25167, 0.0, 1.1)
CellsMeta["Gene_IDs"] <- randomnumbers
head(CellsMeta)
cellsMetaTrim <- subset(CellsMeta, select = c("Gene_IDs"))
integrated_object <-AddMetaData(integrated_object, cellsMetaTrim)
head(integrated_object)

VlnPlot(integrated_object, c("nCount_RNA", "Gene_IDs"), ncol = 2)

DimPlot(integrated_object, group.by = c("orig.ident", "seurat_clusters"))


#integrated_object[["RNA"]] <- split( integrated_object[["RNA"]] ,f =  integrated_object$seurat_clusters )
#################

#integrated_object <- NormalizeData(integrated_object)
#integrated_object <- FindVariableFeatures(integrated_object) 
#integrated_object <-ScaleData(integrated_object)
#integrated_object <-RunPCA(integrated_object)
#integrated_object<-FindNeighbors(integrated_object, dims = 1:30)
#integrated_object <-FindClusters(integrated_object, resolution = 0.05)

#integrated_object <-RunUMAP(integrated_object, dims = 1:30)


#meta.all <- integrated_object %>% 
#  dplyr::rename(Integrated_object = "integrated_object")  %>%
#  select(ID_prefix, sample_ID, Cust_ID, Exp_batch)

##############################################

integrated.clusters <- AddModuleScore(object = integrated.strain, features = microglia.gene.list, ctrl = 100, name = 'Microgial Features')

object.significance.PCA.scores <- JackStraw(integrated.clusters, num.replicate = 30, dims = 30)
object.significance.PCA.scores <-ScoreJackStraw(object.significance.PCA.scores,dims = 1:30)


##### Find cluster markers 

integrated_markers <-FindAllMarkers(integrated.strain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 300 )
integrated_markers <- integrated_markers %>% rownames_to_column(var = "symbol")


## save cell metadata and marker info into rda
meta <- integrated.strain@meta.data %>% select(-starts_with("percent.r"))

meta_integrated_markers <-merge(meta, y= integrated_markers, add.cell.idec= c("meta", "marker"), Project = "meta_markers"  )
