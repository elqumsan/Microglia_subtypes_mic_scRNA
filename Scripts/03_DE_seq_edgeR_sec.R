
library(tidyverse)
library(cowplot)
library(Seurat)
library(extrafont)
library(stringr)
library('tidyr')
library(edgeR)
library(gridExtra)
library("scater")
library(org.Mm.eg.db) ### add gene names

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

orig.symbol  <- rownames(singleCell_object)
#rownames(singleCell_object) == EnsembleID$symbol

gridExtra::grid.arrange(plotUMAP(singleCell_object , color_by = "strain", text_by = "seurat_clusters"))

ggsave(paste(global_var$global$path_DE_seq_edgeR,  "SingleCellExper.png", sep ="/"), width = 3.5, height = 5.3 , units = "in" , dpi = 600 )


### using 'label' and 'sample' as our two factors; each column of the output 
### corresponds to one unique combination of these two factors 
summed <- aggregateAcrossCells(singleCell_object,
                     id= DataFrame(
                                    Lable= singleCell_object$seurat_clusters,
                                   sample = colnames(singleCell_object)))
#label = c("0","1", "2","3","4", "5", "6", "7", "8")
#current <- summed[, label == summed$Lable]
### Below are the testing to make DE gene list for making comparison between the two strains we have 

## Creating up a DGEList object for use in edgeR 

y <- DGEList(counts(summed), samples= colData(summed))

### Add gene names
y$genes$orig_symbol <- orig.symbol

#y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y), keytype = "ENSEMBL"  ,column = "ENTREZID")

#### remove samples with low library size 
discarded <- isOutlier(y$samples$lib.size, log = TRUE, type = "lower")
y <-y[, !discarded]
summary(discarded)


### remove genes that are lowly expressed
keep <- filterByExpr(y, group = summed$strain) ## check group argument when filtering
y <- y$genes[keep ]
summary(keep)

### Trimmed means of M-values of methods for normalization 
y <-calcNormFactors(y)

### Statistical Modeling
y$samples$strain = factor(y$samples$strain, levels = c("Veh", "AZT"))

str(y$samples)


design <- model.matrix(~strain, y$samples)

#y <- estimateDisp(y, design, trend.method = "none", subset = 5000 )

y <- estimateDisp(y )


summary(y$trended.dispersion)

png(filename = paste( global_var$global$path_DE_seq_edgeR,  "PlotBCV.png", sep= "/"))
p <- plotBCV(y)
dev.off()


# ggsave(filename =  paste(global_var$global$path_DE_seq_edgeR, "estimateDisp.png", sep = "/"), print(p),width = 10.5, height = 7)
# ggsave(plot = p,  paste(global_var$global$path_DE_seq_edgeR, "estimateDisp.png", sep = "/"), width = 10.5, height = 7 )


fit <- glmQLFit(y, design , robust = TRUE)

summary(fit$var.prior)
summary(fit$df.prior)

png(filename = paste(global_var$global$path_DE_seq_edgeR, "PlotQLDisp.png", sep= "/") )
p <- plotQLDisp(fit)
dev.off()

#ggsave( plot = p ,paste(global_var$global$path_DE_seq_edgeR, "glmQLFit.png", sep = "/"), width = 3.5, height = 5, units = "in", dpi = 300 )
#dev.off()


colnames(coef(fit))
res <- glmQLFTest(fit, coef = 1)
summary(decideTests(res))

summary(decideTests(res))

x <- topTags(res, n=NULL, p.value= 1)


