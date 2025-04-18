
library(tidyverse)
library(cowplot)
library(Seurat)
library(extrafont)
library(stringr)
library('tidyr')
library(data.table)

# 1. QC for merged data

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

####### Define color palettes and plot themes

colLib = c("blue", "darkorange")   # color for libraries
names(colLib) <- c("Veh", "AZT")

colGEX = c("grey85", brewer.pal(2, "Reds"))
colCcy = c("black", "blue", "darkorange")
plotTheme <- theme_classic(base_size = 18)

inpUMIs <- NULL

####
############# Perform QC - Compute cell QC metrics

#Vehdata.path <- "/shared/home/mabuelqumsan/rnaseqmva_shared_space/TANG_Lab/Xin_data/Veh/raw_feature_bc_matrix.h5"

#AZTdata.path <-   "/shared/home/mabuelqumsan/rnaseqmva_shared_space/TANG_Lab/Xin_data/AZT/raw_feature_bc_matrix.h5"

#Vehdata  <- Read10X_h5( Vehdata.path)
#AZTdata <- Read10X_h5( AZTdata.path)

#inpUMIs <- cbind(Vehdata, AZTdata)


oupQCcell <- data.table(
            sampleID = colnames(integrated.strain@assays$RNA$scale.data),
            library = tstrsplit(colnames(integrated.strain@assays$RNA$scale.data),"_")[[1]],
            nUMI = colSums(integrated.strain@assays$RNA$counts),
            nGene = rowSums(integrated.strain@assays$RNA$counts ), 
            pctMT = colSums(integrated.strain@assays$RNA$counts[grep("^MT-",rownames(integrated.strain@assays$RNA)),] )
)

oupQCcell$pctMT <- 100 * oupQCcell$pctMT / oupQCcell$nUMI
oupQCcell$library <-  factor(oupQCcell$library, levels = names(colLib))

########



### Plot cell QC metrics
p1 <- ggplot(oupQCcell, aes(nUMI, fill = library)) + 
  geom_histogram(binwidth = 500, color = "black") + xlim(c(0,40e3)) + 
  geom_vline(xintercept = c(1.5e3,30e3), color = "red", linetype = "dashed")+ 
  xlab("No. UMIs") + scale_fill_manual(values = colLib) + plotTheme

p2 <- ggplot(oupQCcell, aes(nGene, fill = library)) + 
  geom_histogram(binwidth = 100, color = "black") + xlim(c(0,6e3)) + 
  geom_vline(xintercept = c(0.5e3), color = "red", linetype = "dashed")+ 
  xlab("No. Detected Genes") + scale_fill_manual(values = colLib) + plotTheme

p3 <- ggplot(oupQCcell, aes(pctMT, fill = library)) + 
  geom_histogram(binwidth = 0.5, color = "black") + xlim(c(0,20)) + 
  geom_vline(xintercept = c(15), color = "red", linetype = "dashed")+ 
  xlab("Percentage MT Genes") + scale_fill_manual(values = colLib) + plotTheme


ggsave(p1 + p2 + p3 + guide_area() + plot_layout(nrow = 2, guides = "collect"), 
       width = 10, height = 8, filename = "../Microglia_subtypes_mic_scRNA/findings/01a_QC_strains/basicCellQC.png")

###################

dim(meta)
head(meta)

meta_tidy <- meta %>% 
  select(orig.ident,nCount_RNA, nFeature_RNA, Strain) %>%
  gather(-orig.ident, -Strain ,key= "QC", value="value" )

head(meta_tidy)

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
