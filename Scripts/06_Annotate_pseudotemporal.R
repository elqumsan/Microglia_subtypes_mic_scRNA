install.packages("devtools")
devtools::install_github("immunogenomics/presto", force = TRUE)
library(presto)
library(dplyr)

devtools::install_github('quadbiolab/voxhunt')
library(voxhunt)

####### Differrential expression analaysis with DESeq2


counts <- integrated.strain@assays$RNA$counts
metadata <- integrated.strain@meta.data

metadata$cluster_id <- factor(integrated.strain@active.ident)

sce <- SingleCellExperiment(assays = list(counts = counts),
                     colData = metadata)
assay(sce)
dim(counts(sce))
counts(sce) [1:6, 1:6]
dim(colData(sce))
head(colData(sce))

cluster_names <- levels(colData(sce)$cluster_id)
length(cluster_names)

groups <- colData(sce)[, c("cluster_id", "strain")]
head(groups)


aggr__counts <- aggregate(t(counts(sce)), groupings = groups, fun = "sum" )


######### Annotate cell clusters analysis


cl_markers_presto <-wilcoxauc(integrated.strain)

cl_markers_presto %>%  
  filter( logFC > log(1.2) & pct_in > 20 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group = T) %>%
  top_n(n=2, wt = logFC) %>%
  print(n= 40, width = Inf)

object.ref <- IntegrateLayers(object = integrated_object, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose= FALSE)
object.ref<- FindNeighbors(object.ref, reduction = "integrated.cca", dims = 1:30)
object.ref <- FindClusters(object.ref)
RunUMAP(object.ref, reduction = "integrated.cca", dims= 1:30)
DimPlot(object.ref )

######### visualizing the identified top cluster markers 
cl_markers <- FindAllMarkers(integrated.strain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC)

DoHeatmap(integrated.strain, features = top10_cl_markers$gene) + NoLegend()

plot1 <- FeaturePlot( integrated.strain,  c("Cst7",  "Tmem119"), order = T)
plot2 <- VlnPlot(integrated.strain, features =  c("Cst7","Tmem119", "Ifit3"), pt.size = 0)




####################
variable_genes('Tmem119', 200)$gene

#################################

genes <- c("Cst7",  "Apoe", "Cx3cr1", "Tmem119", "Ifit3", "Ifitm3", "Irf7", "Hexb", "Cd81", "Cst3", "Rplp1", "Rps21", "Rps24" 
           , "C3ar1", "Stmn1", "Top2a", "Birc5" )

DoHeatmap(integrated.strain, features = genes) + NoLegend()

integrated_markers %>% group_by(cluster) %>% top_n(n= 2, wt = avg_log2FC)
