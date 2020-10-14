lapply(c('tidyverse', 'cowplot', 'gridExtra','patchwork', 'scales', 'biomaRt', 'pheatmap'), require, character.only = TRUE)

sample_tbl <- read_tsv('data/obds_sampletable.tsv')
counts_tbl <- read_tsv('data/obds_countstable.tsv.gz')

# listMarts()
# useMart(biomart='ENSEMBL_MART_ENSEMBL') %>% 
#     listDatasets() %>%
#     # View()
#     filter(grepl('^mmus', dataset))
# 
# useMart(biomart='ENSEMBL_MART_ENSEMBL',
#         dataset='mmusculus_gene_ensembl') %>% 
#     listFilters() %>% 
#     View()
# 
# useMart(biomart='ENSEMBL_MART_ENSEMBL',
#         dataset='mmusculus_gene_ensembl') %>% 
#     listAttributes() %>% 
#     View()

conversion_tbl <- useMart(biomart='ENSEMBL_MART_ENSEMBL',
               dataset='mmusculus_gene_ensembl') %>% 
    getBM(mart=., attributes=c('mgi_symbol', 'ensembl_gene_id')) %>%
    as_tibble()

# tidy the counts table and add the gene symbol
counts_tdy <- counts_tbl %>%
    pivot_longer(-Geneid, names_to = 'sampleid', values_to = 'count') %>% 
    left_join(conversion_tbl, by=c('Geneid'='ensembl_gene_id'))

# tidy the sample metadata table
sample_tdy <- sample_tbl %>% 
    select(-c(library_layout, species)) %>% 
    separate(sample_title, c('gene', 'genotype', 'celltype', 'replicate'), sep='_') %>% 
    unite('Genotype', gene, genotype)

# add metadata to counts table
combine_tdy <- left_join(counts_tdy, sample_tdy, by=c('sampleid' = 'Sample_accession'))

# convert to CPM with log conversion
processed_tdy <- combine_tdy %>% 
    group_by(sampleid) %>% 
    mutate(total_count_per_sample = sum(count)) %>% 
    mutate(total_count_in_million = total_count_per_sample / 1e6) %>% 
    mutate(CPM = count / total_count_in_million) %>%
    mutate(log2CPM = log2(CPM + 1)) %>% 
    select(-c(total_count_in_million, total_count_per_sample))

# Plot read depth per sample
processed_tdy %>% 
    summarise(read_depth=sum(count)) %>% 
    ggplot(aes(x=sampleid, y=read_depth)) +
    geom_bar(fill='lightblue', stat='identity') +
    labs(title='Read depth per sample',
         x='Sample ID',
         y='Read Depth') +
    theme(plot.title = element_text(hjust=0.5),
          axis.text.x = element_text(angle=90, vjust=0.5))

# How many genes have no counts for any sample?
processed_tdy %>% 
    group_by(Geneid) %>% 
    summarise(total_counts=sum(count)) %>% 
    filter(total_counts==0) %>% 
    nrow()

# Draw a density plot of log2(CPM + 1) for all genes
processed_tdy %>% 
    ggplot(aes(x=log2CPM)) +
    geom_density(aes(color=sampleid))
    
# Filter out genes that have low expression in 3 or more samples (ie CPM < 0.5)
filtered_tdy <-  processed_tdy %>% 
    group_by(Geneid) %>% 
    mutate(low_expression=sum(CPM<0.5)) %>% 
    filter(low_expression<3) %>% 
    ungroup

# What proportion of genes are lowly expressed?
genes_total <- length(unique(processed_tdy$Geneid))
genes_remain <- length(unique(filtered_tdy$Geneid))
(genes_total - genes_remain) / genes_total

# Make a density plot of log2(CPM + 1) with the filtered data
filtered_tdy %>% 
    ggplot(aes(x=log2CPM)) +
    geom_density(aes(color=sampleid))

# Plot CD4 and CD8 expression for all samples, colour by replicate and facet by genotype against cell type
filtered_tdy %>% 
    filter(mgi_symbol == 'Cd4' | mgi_symbol == 'Cd8a' | mgi_symbol == 'Cd8b1') %>% 
    ggplot(aes(x=log2CPM, y=mgi_symbol, color=replicate)) +
    geom_boxplot() +
    facet_grid(celltype ~ Genotype)

# Generate the same plot for Egr2 and Egr3 for all samples
filtered_tdy %>% 
    filter(mgi_symbol == 'Egr2' | mgi_symbol == 'Egr3') %>% 
    ggplot(aes(x=log2CPM, y=mgi_symbol, color=replicate)) +
    geom_boxplot() +
    facet_grid(celltype ~ Genotype)

# Choose 8 biologically relevant genes and plot a heatmap using the pheatmap package
relevant_genes <- c('Dnmt1','Dnmt3a','Dnmt3b','Tet1','Tet2','Tet3','Arid1a','Arid1b')
filtered_tdy %>% 
    filter(mgi_symbol %in% relevant_genes) %>%
    select(c(mgi_symbol, Genotype, celltype, replicate, log2CPM)) %>%
    pivot_wider(names_from=c(Genotype, celltype, replicate), names_sep='_', values_from=log2CPM) %>% 
    arrange(mgi_symbol) %>% 
    column_to_rownames('mgi_symbol') %>%
    pheatmap
