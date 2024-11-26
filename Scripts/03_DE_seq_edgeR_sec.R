
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


### Single cell RNAseq is notorious for having a normal distribution for gene expression
# scRNAseq has a trasnscript recovery rate between 60%-80%

png(filename = paste( global_var$global$path_DE_seq_edgeR,  "density_of_Graph_2500_randome_cells.png", sep ="/"))

plot(density(sample(JoinLayers(integrated_object@assays$RNA)$count["Gapdh",],2500)),cex=0,lty=1, main="Density of Gapdh in 2500 random cells")
dev.off()

png(filename = paste( global_var$global$path_DE_seq_edgeR,  "Histogram_of_graph.png", sep ="/"))
hist(sample(JoinLayers(integrated_object@assays$RNA)$count["Gapdh",],2500),breaks=99,main="Histogram of Gapdh in 2500 random cells",ylab="Frequency",xlab="Gene counts")
dev.off()


#########
dim(integrated_markers)
integrated_markers <- integrated_markers[integrated_markers$p_val_adj < 0.05,]
integrated_markers <- integrated_markers[order(integrated_markers$avg_log2FC, decreasing = T ), ]

png(filename = paste(global_var$global$path_DE_seq_edgeR, "Heatmap.png", sep = "/"))

DoHeatmap(object = integrated.strain, features = (integrated_markers$gene[1:50]), group.by = "strain", slot = "data", assay = "RNA")
dev.off()

png(filename = paste(global_var$global$path_DE_seq_edgeR, "DimHeatmap.png" , sep = "/"))
DimHeatmap(integrated.strain, dims = 1:9, cells = 500 , balanced = TRUE, ncol = 3)
dev.off()

png(filename = paste(global_var$global$path_DE_seq_edgeR, "FeaturePlot.png", sep = "/" ))
FeaturePlot(object = integrated.strain, features = head( integrated_markers$gene ))
dev.off()



##### Feature selection for following heterogeneity analysis
feature_variable <- FindVariableFeatures(integrated.strain, nfeatures = 3000)
top_feature <- head(VariableFeatures(feature_variable), 20)

plot1 <- VariableFeaturePlot(feature_variable)
plot2 <- LabelPoints(plot = plot1 , points = top_feature, repel = TRUE)
x <- plot1 + plot2
ggsave(filename = paste(global_var$global$path_DE_seq_edgeR, "heterogeneity_analysis.png", sep = "/"), units = "in" , width = 10, height = 5, dpi = 300 )



########## Non-linear dimension reduction for visualization 
object1 <- RunTSNE(integrated.strain, dims = 1:20)
object2 <- RunUMAP(integrated.strain, dims = 1:20)

plot1 <- TSNEPlot(object1)
plot2 <-UMAPPlot(object2)

plot1 + plot2

ggsave(filename = paste(global_var$global$path_DE_seq_edgeR, "tSNE_UMAP_plot.png", sep = "/"), units = "in", width = 10, height = 5, dpi = 300 )

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

write_delim(x$table, paste(global_var$global$path_DE_seq_edgeR, "glmQLFTest.txt", sep="/"), delim = "\t")
