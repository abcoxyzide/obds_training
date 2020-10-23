sapply(c('tidyverse', 'DESeq2', 'cowplot', 'EnsDb.Mmusculus.v79', 'vsn', 'pheatmap', 'ggrepel'), require, character.only = T)
if (!requireNamespace('patchwork')) {install.packages('patchwork'); require('patchwork', character.only = T)}

df_sample <- read_tsv('data/obds_sampletable.tsv')
df_count <- read_tsv('data/obds_countstable.tsv.gz')

data.m <- df_count %>% column_to_rownames('Geneid') %>% as.matrix
data.sample <- df_sample %>% 
    separate(sample_title, c('locus', 'genotype', 'celltype', 'replicate'), sep='_') %>% 
    unite(col='condition', locus, genotype, celltype, sep='_') %>% 
    dplyr::select(-c(species, library_layout)) %>% 
    mutate(condition = factor(condition, levels=c("Egr2/3_DKO_CD8", "Egr2/3_DKO_CD4", "Egr2_Kin_CD8", "Egr2_Kin_CD4"))) %>% 
    column_to_rownames('Sample_accession')

# neat way of checking how many T/F came out from that statement
table(colnames(data.m)==rownames(data.sample))

# making dds object
dds <- DESeqDataSetFromMatrix(data.m, data.sample, ~condition)

colData(dds) # or dds@colData
rowRanges(dds)
design(dds)
counts(dds) # equivalent to assays(dds)$counts

#-------------
# # estimate the size factors, and return the object to itself
# dds <- estimateSizeFactors(dds)
# # estimate the dispersion
# dds <- estimateDispersions(dds)
# # perform wald test
# dds <- nbinomWaldTest(dds)
#-------------
# a wrapper for the above 3 functions
dds <- DESeq(dds)

# see the size factor
sizeFactors(dds) 
sizeFactors(dds) %>% 
    data.frame(size=. ,
               group=colData(dds)$condition) %>%
    ggplot(aes(x=rownames(.), y=size)) +
    geom_col(aes(fill=group)) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))

# see the dispersion
dispersions(dds) 
plotDispEsts(dds) # looks good cus similar to example they've given??

# see the coefficients after Wald test
head(coef(dds))

#-------------
# DESeq saves the result to a specific data format
res <- results(dds, contrast=c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'))
res <- results(dds, name = 'condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8')


p_padj <- 
    res %>% 
    as.data.frame %>% 
    ggplot(aes(x=padj))+
    geom_histogram(color='black', fill='yellow', alpha=0.5, )
p_pv <- 
    res %>% 
    as.data.frame %>% 
    ggplot(aes(x=pvalue))+
    geom_histogram(color='black', fill='blue', alpha=0.5)

# p_padj / p_pv
# 
# try(require(patchwork))


# compare CD4 vs CD8 DKO (as specified above with the 'contrast' argument)
plotMA(res)


# compare different methods of shrinkage to visualize stuff
# coef is essentially the same as contrast; but apeglm annoyingly needs the coef option
resNorm <- lfcShrink(dds, contrast=c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type ='normal')
resAsh <- lfcShrink(dds, contrast=c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type ='ashr')
resApeg <- lfcShrink(dds, coef='condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8', type ='apeglm') # default

# compare using MA plots
# but it's a base R function, and cannot be saved directly to an object, thus needs recordPlot() if you want to use function like plot_grid
# alternatively, just use par(mfrow=c(n, m))
plotMA(resNorm) ; resNormplot <- recordPlot()
plotMA(resAsh) ; resAshplot <- recordPlot()
plotMA(resApeg) ; resApegplot <- recordPlot()

plot_grid(resNormplot, resApegplot, resApegplot)

par(mfrow=c(2,2)); plotMA(res); plotMA(resNorm); plotMA(resAsh); plotMA(resApeg)

# using the EnsDb.Mmusculus.v79 package to get gene names from gene id
# manual: https://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
edb <- EnsDb.Mmusculus.v79
columns(edb)
keytypes(edb)

resApeg_df <- as.data.frame(resApeg)

conv_df <- select(edb, keys= rownames(resApeg_df), columns = c('GENEID', 'SYMBOL'), keytype = 'GENEID')

resApeg_df <- resApeg_df %>% 
    rownames_to_column('GENEID') %>% 
    left_join(conv_df)

# remove those entries with no mappable genes
resApeg_padj <- resApeg_df %>% .[!is.na(.$SYMBOL), ] %>% .[!is.na(.$padj), ]
write_csv(resApeg_padj, 'data/results_apeg.csv')


# get significant genes and log count > 1
sig_res <- resApeg_padj %>% .[.$padj <0.05, ] %>% .[.$log2FoldChange > 1, ] %>% arrange(-abs(log2FoldChange))
top20 <- head(sig_res, 20)

# calculate vst and plot meansdplot
vst <- vst(dds, blind=F)
rlog <- rlog(dds, blind=F)
assay(vst)

meanSdPlot(assay(vst))
meanSdPlot(assay(rlog))

# plot pca
plotPCA(vst, ntop=nrow(vst)) # doing good
plotPCA(vst, ntop=500) # doing even better; because reduced noise


# heatmap of top 20 genes
vst_top20 <- assay(vst) %>% .[rownames(.) %in% top20$GENEID, ]
vst_top20 %>% pheatmap(scale='row')
    

# volcano plot
# y = -log10(pvalue), x = fold change
df_for_volcanoplot <- resApeg_padj %>% 
    dplyr::filter(!is.na(padj))


color <- c('black', 'red')
highlight <- with(df_for_volcanoplot, factor(padj<0.05 & abs(log2FoldChange)>1))
color <- color[highlight]

df_for_volcanoplot %>% 
    mutate(y=-log10(padj)) %>% 
    mutate(highlight=(padj<0.05 & abs(log2FoldChange)>1)) %>% 
    ggplot(aes(x=log2FoldChange, y=y)) +
    geom_point(color=color)

genes_of_interest <- c('Atp8b4', 'Ccdc7', 'Fscn1', 'Kcnj8', 'Rasd2', 'Sox8')

df_for_volcanoplot <- df_for_volcanoplot %>% 
    mutate(y=-log10(padj)) %>% 
    mutate(highlight=(padj<0.05 & abs(log2FoldChange)>1))
ggplot(df_for_volcanoplot, aes(x=log2FoldChange, y=y, label=SYMBOL)) +
    geom_point(color=color) +
    geom_label_repel(data = df_for_volcanoplot[df_for_volcanoplot$SYMBOL %in% genes_of_interest, ])



