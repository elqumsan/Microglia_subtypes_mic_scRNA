
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
emapplot(pairwise.prop.test(ora_results))
emapplot()
