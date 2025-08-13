library(biomaRt)


##### Prepare marker gene list 
integrated_object <- JoinLayers(integrated_object)
markers  <- FindAllMarkers(integrated_object, only.pos = TRUE, min.pct = 0.25)
filtered_markers <- markers %>% 
  filter(p_val_adj < 0.05 , avg_log2FC > 0.25 )

##### Obtain MSigDB (Mouse) gene sets
### Identify mouse-specific gene sets related to microglia from MSigDB . use the C7 (immune signatures), C8 (cell type signatures)
### Get Microglia-related cell type gene sets from mouse 
msigdb_mouse <- msigdbr(species = "Mus musculus", category = "C8")
microglia_sets <- msigdb_mouse %>% 
                  filter(grepl("Microglia", gs_name, ignore.case = TRUE))


msig_df <- microglia_sets %>% 
          dplyr::select(gs_name, gene_symbol) %>%
          as.data.frame()

##### Performing over-representaation Analysis (ORA)
### ORA tests whetther your marke gene list is enriched in known Microglia gene sets (from MSigDB) 

### Extract your gene list
gene_list <- unique(filtered_markers$gene)


### Convert to SYMBOL
entrez_genes <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

### Prepare MSigDB list
# msig_list <- split(microglia_sets$gene_symbol, microglia_sets$gs_name)


##### Performing Over-Representation Analysis (ORA)
### Run ORA
ora_results <- enricher(entrez_genes$SYMBOL, TERM2GENE = msig_df)

head(ora_results)

## Dotplot of enriched terms
dotplot(ora_results, showCategory = 20, font.size = 12)

## Barplot 
barplot(ora_results, showCategory = 15 , title = "Microglai Signature Enrichment")

### Enrichment map (if multiple terms enriched)
## emapplot(pairwise.prop.test(ora_results))
## emapplot()

####### Mouse cell cycle Genes 
######### 
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
## Load Human Seurat cell cycle genes
s_genes <- cc.genes.updated.2019$s.genes
g2m_genes <- cc.genes.updated.2019$g2m.genes

#### convert to mouse genes names
convert_cc_genes  <- function(genes){
  getLDS(attributes = c("hgnc_symbol") , filters = "hgnc_symbol", 
         values = genes, mart = human,
         attributesL = c("mgi_symbol"), martL = mouse,
         uniqueRows = TRUE)$MGI.symbol 
}

s_genes_mouse <- convert_cc_genes(s_genes)

integrated_object <- CellCycleScoring(integrated_object, s.features =  s_genes,
                 g2m.features = g2m_genes ,
                 set.ident = FALSE)


DimPlot(integrated_object, group.by = "Phase", label = TRUE)

VInPlot