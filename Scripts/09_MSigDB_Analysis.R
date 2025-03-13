library(msigdb)
library(ExperimentHub)
library(GSEABase)

## Download data from the msigdb 

eh <-ExperimentHub(ask = FALSE) 
query(eh, 'msigdb')

eh[['EH5421']]

##### Accessing the Human MsigDB
msigdb.hs <- msigdb::getMsigdb(org = 'hs', id = 'SYM' )

### Accessing the mouse MSigDB
msigdb.mm <- msigdb::getMsigdb(org = "mm", id = 'SYM')

#### Download and integrating KEGG gene sets 
msigdb.hs= appendKEGG(msigdb.hs)

msigdb.mm = appendKEGG(msigdb.mm)

msigdb.hs
msigdb.mm

### Acessing the GeneSet and genesetCollection objects 
length(msigdb.hs)

length(msigdb.mm)


msigdb.mm[[1000]]
