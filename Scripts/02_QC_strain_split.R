
library(tidyverse)
library(cowplot)
library(Seurat)
library(extrafont)
library(stringr)
library('tidyr')
library(data.table)
library(dplyr)

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

integrated.strain <- integrated_object
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
######  Remove Low quality cells ()

oupQCcell <- oupQCcell[nUMI > 1.5e3]
oupQCcell <- oupQCcell[nUMI < 30e3]
oupQCcell <- oupQCcell[nGene > 0.5e3]
oupQCcell <- oupQCcell[pctMT < 15]

integrated.strain <- JoinLayers(integrated.strain) 


inpUMIs <- integrated.strain@assays$RNA$counts[,as.character(oupQCcell$sampleID)]

###########
### Compute gene QC metrics
inpUMIs <- inpUMIs[rowSums(inpUMIs) > 0, ]  # Remove genes with no reads

oupQCgene <- data.table(
              gene = rownames(inpUMIs),
              meanRead = rowSums(inpUMIs) / ncol(inpUMIs),
              cellExpr = rowSums(inpUMIs > 0),
              library = tstrsplit(colnames(inpUMIs),"_")[[1]]
)

oupQCgene$log10meanRead <- log10(oupQCgene$meanRead)
oupQCgene$log10cellExpr <- log10(oupQCgene$cellExpr)

### Plot gene QC metrics
p1 <- ggplot(oupQCgene, aes(log10meanRead, fill = library)) + 
  geom_histogram(binwidth = 0.05, color = "black") +
  xlab("log10(Average UMIs)") + plotTheme

p2 <- ggplot(oupQCgene[cellExpr != 0], aes(log10cellExpr, fill = library) ) + 
  geom_histogram(binwidth = 0.05 , color = "black") + 
  geom_vline(xintercept = log10(8), color = "red", linetype = "dashed") +
  xlab ("log10(No. Cells Expression gene)" ) + plotTheme

ggsave(p1 +p2, width = 10 , height = 4, filename = "../Microglia_subtypes_mic_scRNA/findings/01a_QC_strains/basicGeneQC.png")

#### Remove Lowly expressed genes (23130 genes down to 19064 genes post - QC )

oupQCgene <- oupQCgene [cellExpr >= 8]
inpUMIs <- inpUMIs[as.character(oupQCgene$gene),]

########## Create Seurat object + Preprocessing 
### Create Seurat Object
### nCount_RNA = & no. UMI count ; nFeature_RNA = no. detected genes
### We will also compute percentage MT genes here and set the library

seu <- CreateSeuratObject(inpUMIs, project = "Just Microglia")
colnames(seu@meta.data)
seu$pctMT <- 100 * colSums(inpUMIs[grep("^MT-", rownames(inpUMIs)),])

seu$pctMT <- seu$pctMT / seu$nCount_RNA

#### Create library and donor column
seu$library = factor(seu$orig.ident, levels = names(colLib))
Idents(seu) <- seu$library
seu$donor <- gsub("pos|neg" , "" , seu$library)
seu$donor = factor(seu$donor)

##### Normalization + Feature selection
seu <- NormalizeData(seu, assay = "RNA")
seu <-FindVariableFeatures(seu, selection.method = "vst"  )

top10 <- head(VariableFeatures(seu), 20)
##### Plot gene expression mean vs variance
ggData = as.data.frame( as.matrix( seu[["RNA"]]@meta.data))

# ggData = as.data.frame(seu@assays$RNA)

p1 <- VariableFeaturePlot(seu)
p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE) 
p1 + p2 

ggsave(p1 + p2 + plot_layout(guides = "collect"),
       width = 20, height = 10, filename =  "../Microglia_subtypes_mic_scRNA/findings/01a_QC_strains/basicHVG.png")

#p1 <- ggplot(ggData, aes(vf_vst_counts_mean, vf_vst_counts_variance, color = vf_vst_counts_variable)) + 
#  geom_point() + scale_color_manual(values = c("black", "firebrick")) + 
#  geom_line(aes(vf_vst_counts_mean, vf_vst_counts_variance.expected), color = "dodgerblue") + 
#  xlab("Average Expression") + ylab("Gene Variance") +
#    plotTheme
 
####### PCA / tSNE / UMAP dimension reduction
seu <- ScaleData(object = seu)  ### Scale data prior to PCA
seu <- RunPCA(seu)
 
p1 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE,
               cols = colLib) + plotTheme + coord_fixed()  
p2 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE,
              cols = colLib, dims = c(1,3)) + plotTheme + coord_fixed()


###### Elbow plot
p1 <- ElbowPlot(seu, ndims = 40)
ggData <- seu@reductions[["pca"]]@stdev  # fit lines for elbo plot
lm(y~x, data = data.frame(x = 16:20 , y = ggData[16:20]))
lm(y~x , data = data.frame(x= 31:35 , y=ggData[31:35]) )
nPC <- 20 ## determined from elbow plot
p1 <- p1 + plotTheme + geom_vline(xintercept = nPC, color = "firebrick") + 
  geom_path(data = data.frame(x= 16:35, y = 16:35*-0.07879 + 3.84064) ,
            aes(x,y) , color = "dodgerblue") + 
  geom_path(data = data.frame(x = 16:35 , y= 16:35* -0.02506 + 2.64554) , 
            aes(x,y) , color = "dodgerblue")

ggsave(p1 , width = 5 , height = 4 , filename = "../Microglia_subtypes_mic_scRNA/findings/02_QC_strain_split/basicPCAelbow.png")

## Run tSNE 
seu <- RunTSNE(seu , dims = 1:nPC , num_threads = nPC)
p1 <- DimPlot(seu, reduction = "tsne", pt.size = 0.1, shuffle = TRUE, cols = colLib) + plotTheme + coord_fixed()


### RUN UUMAP
seu <- RunUMAP(seu, dims = 1:nPC)
p2 <- DimPlot(seu ,  reduction = "umap", pt.size = 0.1 , shuffle = TRUE,
              cols =  colLib) + plotTheme + coord_fixed()

ggsave(p1 + p2 + plot_layout(guides = "collect"), width = 10 , height = 4 , filename = "../Microglia_subtypes_mic_scRNA/findings/02_QC_strain_split/basicUmapGex.png")

### plot some genes from literature
p1 <- FeaturePlot(seu, reduction = "pca", pt.size = 0.1,
            features =  c( "Fcrls","P2ry12","Cx3cr1", "Trem2", "C1qa" ,"Tmem119", "Tubb3", "Meg3", "Dcx", "Gfap", "Olig1", "Vtn", 
                           "F13a1"), order = TRUE) &
  scale_color_gradientn(colours = colGEX) & plotTheme & coord_fixed()
#############
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
