require(GO.db)
require(EnsDb.Mmusculus.v79)
require(clusterProfiler)
require(tidyverse)

# documentation for the package
# http://yulab-smu.top/clusterProfiler-book/chapter12.html


data <- read.csv('data/results_apeg.csv')

# map entrez id
keys <- data$GENEID
# columns(EnsDb.Mmusculus.v79)
# keytypes(EnsDb.Mmusculus.v79)
geneid2andrezid <- ensembldb::select(EnsDb.Mmusculus.v79,
                      keys=keys,
                      columns=c("GENEID", "ENTREZID"),
                      keytype="GENEID")
data <- data %>% 
    left_join(geneid2andrezid)

# upregulated = padj < 0.05 & log fold change > 1
# downregulated = log fold change < -1
up_genes <- data %>% 
    filter(!is.na(ENTREZID)) %>% 
    filter(padj < 0.05) %>% 
    filter(log2FoldChange > 1)

dn_genes <- data %>% 
    filter(!is.na(ENTREZID)) %>% 
    filter(padj < 0.05) %>% 
    filter(log2FoldChange < -1)

# Go db requires to supply a named vector of log fold changes
genelist_up <- with(up_genes, set_names(log2FoldChange, ENTREZID)) %>% sort(decreasing = T)
genelist_dn <- with(dn_genes, set_names(log2FoldChange, ENTREZID)) %>% sort(decreasing = F)
genelist_all <- data %>% 
    filter(!is.na(ENTREZID)) %>% 
    pull(log2FoldChange, name=ENTREZID)

#-----    
## perform ORA

# with GO
ora_go <- enrichGO(gene = names(genelist_up),
                universe = names(genelist_all),
                OrgDb = 'org.Mm.eg.db')

# with KEGG
# search_kegg_organism('mouse', by='common_name')
# kegg_code scientific_name common_name
#       mmu    Mus musculus       mouse

ora_kegg <- enrichKEGG(names(genelist_up),
                       organism = 'mmu',
                       universe = names(genelist_all))
ora_kegg@result %>% View
ora_kegg_res <- ora_kegg@result


#-----
## visualize

# barplot
# ?enrichplot::barplot.enrichResult
barplot(ora_kegg, showCategory=20) # output 6 pathways; because only 6 are significant
barplot(ora_go, showCategory = 20)

# dotplot
dotplot(ora_go, showCategory = 20)

# Gene-Concept Network
res_symbol <- 
    ora_go %>% 
    setReadable('org.Mm.eg.db', 'ENTREZID')
cnetplot(res_symbol, layout = 'kk')

#-----
## perform GSEA

# with GO
genelist_all_dedup <- sort(genelist_all[!duplicated(names(genelist_all))], decreasing = T)
gse_go <- gseGO(genelist_all_dedup,
                OrgDb = org.Mm.eg.db,
                ont          = "ALL", # BP = biological process, CC = cellular component, MF = molecular function 
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = T
                )
gse_go@result %>% View

## visualization

# ridgeplot
ridgeplot(gse_go)

# gseaplot
p1 <- gseaplot(gse_go, geneSetID = 1, by = "runningScore", title = gse_go$Description[1])
p2 <- gseaplot(gse_go, geneSetID = 1, by = "preranked", title = gse_go$Description[1])
p3 <- gseaplot(gse_go, geneSetID = 1, title = gse_go$Description[1])
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])



