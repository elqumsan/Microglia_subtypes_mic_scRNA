
library(tidyverse)
library(cowplot)
library(Seurat)
library(extrafont)
library(stringr)
library('tidyr')
library(edgeR)
library(gridExtra)
library("scater")

# convert gene names to gene ID in the row name of integrated object
## This script is to determine DE gene affected by strain
## Use EdgeR QFT method adapted from : Tutorial https://osca.bioconductor.org/multi-sample-comparisons.html#motivation-8
# Bench mark DE gene analysis identified EdgeR QFT as one of the top methods  https://www.nature.com/articles/nmeth.4612


#out_path <- "../Microglia_subtypes_mic_scRNA-/findings/02_QC_strain_split/"

#Wtdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/Veh/")
#WT_Project <- "WT_data_Microglia"
#AZTdata.path <- ("/shared/ifbstor1/projects/rnaseqmva/TANG_Lab/Xin_data/AZT/")
#AZT_Project <-  "AZT_data_Microglia"

#WTdata  <- Read10X(data.dir = Wtdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
#AZTdata <- Read10X(data.dir = AZTdata.path, gene.column = 2, cell.column = 1, unique.features = T, strip.suffix = F)
#WT_object <-  CreateSeuratObject(counts = WTdata, project = "Microglia subtypes_WT", min.cells = 3, min.features = 300)
#AZT_object <- CreateSeuratObject(counts = AZTdata, project = "Microglia subtypes_AZT", min.cells = 3, min.features = 300)

#integrated_object <-  merge(WT_object, y = AZT_object, add.cell.ids = c("WT", "AZT"), merge.data= TRUE)

#integrated_object <- integrated_object %>% NormalizeData()  %>% FindVariableFeatures()  %>% ScaleData()  %>% RunPCA() %>% 
#                   FindNeighbors(dims = 1:30 , reduction = "pca") %>%
#                   FindClusters(resolution = 0.05) 
  


  
              
# integrated_object <-integrated_object %>% RunUMAP(dims = 1:30 , reduction.name = "umap")



#integrated_object <- RunPCA(integrated_object)
#integrated_object <-FindNeighbors(integrated_object, dims = 1:30)
#cluster_object <- FindClusters(integrated_object,resolution = 0.05)

#strain <- c("Microglia subtypes_WT", "Microglia subtypes_AZT")
#cols =  c("#888888", "#00AA00")

#integrated_object[["Strain"]] <- factor(integrated_object@meta.data$orig.ident, levels = strain)
#integrated_object$Strain <- str_replace(integrated_object$Strain,pattern = "Microglia subtypes_WT", replacement = "WT" )
#integrated_object$Strain <- str_replace(integrated_object$Strain, pattern ="Microglia subtypes_AZT", replacement = "AZT" )
#meta <- integrated_object@meta.data
#dim(meta)
#head(meta)

#meta_tidy <- meta %>% 
#  select(orig.ident,nCount_RNA, nFeature_RNA, Strain) %>%
#  gather(-orig.ident, -Strain ,key= "QC", value="value" )

#head(meta_tidy)


### convert to single cell experiment
DefaultAssay(integrated.strain) <- "RNA"
singleCell_object <- as.SingleCellExperiment(integrated.strain)

#rownames(singleCell_object) == EnsembleID$symbol

gridExtra::grid.arrange(plotUMAP(singleCell_object , color_by = "strain", text_by = "seurat_clusters"))

ggsave(paste(global_var$global$path_DE_seq_edgeR,  "SingleCellExper.png", sep ="/"), width = 3.5, height = 5.3 , units = "in" , dpi = 600 )


### using 'label' and 'sample' as our two factors; each column of the output 
### corresponds to one unique combination of these two factors 
summed <- aggregateAcrossCells(singleCell_object,
                     id= DataFrame(
                                    Lable= singleCell_object$seurat_clusters,
                                   sample = colnames(singleCell_object)))


### Below are the testing to make DE gene list for making comparison between the two strains we have 

## Creating up a DGEList object for use in edgeR 

y <- DGEList(counts(summed), samples= colData(summed))

### Statistical Modeling
y$samples$strain = factor(y$samples$strain, levels = c("Veh", "AZT"))

str(y$samples)
### Trimmed means of M-values of methods for normalization 
y <-calcNormFactors(y)

design <- model.matrix(~strain, y$samples)

y <- estimateDisp(y, design )


summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design , robust = TRUE)

summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)

colnames(coef(fit))
res <- glmQLFTest(fit, coef = 4)
summary(decideTests(res))

x <- topTags(res, n=NULL, p.vvalue= 1)



# plot
meta_tidy %>%
  ggplot(aes(y=value, x=orig.ident, color=Strain)) +
  facet_grid(QC ~ Strain, scales = "free_y") +
  geom_violin() +
  geom_boxplot(width= 0.15, outlier.shape = NA, alpha = 0.7 ) +
#  scale_colour_manual(values = cols)
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12, family = "Arial"),
        axis.title.x = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x =  element_blank()
        )
ggsave(filename = "QC01_merged_integrated.png", path = out_path, width = 3.5, height = 4.5, dpi = 300)


##### wrap into function for the following plots

QC_plot <- function(data){
  
  strain <- c("Microglia subtypes_WT", "Microglia subtypes_AZT")
  cols =  c("#888888", "#00AA00")
  
  
  p <- data %>%
    select(orig.ident,nCount_RNA, nFeature_RNA, Strain) %>%
    gather(-orig.ident, -Strain ,key= "QC", value="value" )
    
    ggplot(aes(y=value, x=orig.ident, color=Strain)) +
    facet_grid(QC ~ Strain, scales = "free_y") +
    geom_violin() +
    geom_boxplot(width= 0.15, outlier.shape = NA, alpha = 0.7 ) +
    #  scale_colour_manual(values = cols)
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          strip.text = element_text(face = "bold", size = 12, family = "Arial"),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x =  element_blank()
    ) 
  
    return(p)
}


##### 
  meta %>%
  group_by(Strain) %>%
  summarise(med_nCount_RNA= median(nCount_RNA),
            med_nFeature_RNA=median(nFeature_RNA),
            N_cell = n())


##### Plot WT/AZT expression levels from single-cell RNA-seq regardless of cluster
integrated_object$Genotype <- factor(integrated_object@meta.data$Strain, levels = c("WT", "AZT"))
integrated_object$clusters <- factor(integrated_object$seurat_clusters, levels = c("1","2","3","4","5","6","7"))
VlnPlot(integrated_object, features = c('P2ry12','Cx3cr1',"Ctss",'Tmem119', "Itgam"), pt.size = 0, split.by = "Genotype")
theme(legend.position = "right",
      axis.title = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      title = element_text(size = 10, family = "Arial")
      )

ggsave(filename = "QC_Itgam_single_cell_cluster.png", path = out_path, width = 4, height = 2, dpi = 300)
