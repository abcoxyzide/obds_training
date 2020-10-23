library(tidyverse)
library(biomaRt)
library(pheatmap)
library(ggplot2)
library(DESeq2)
library(EnsDb.Mmusculus.v79)
library(ggrepel)

#Read in data
sample_table <- read_tsv("data/obds_sampletable.tsv")
count_table <- read_tsv("data/obds_countstable.tsv.gz")
# View(sample_table)
# View(count_table)

#Convert the counts table (obds_countstable.tsv.gz) and the sample information table (obds_sampletable.tsv) into a suitable format for generating a DESeqDataSet object

class(count_table)
class(sample_table) #currently both data frames

count_table <- column_to_rownames(count_table, "Geneid")
count_table <- as.matrix(count_table)

sample_table <- column_to_rownames(sample_table, "Sample_accession")

#Check that the samples correspond to one another in both tables
table(colnames(count_table)==rownames(sample_table))

#Set Egr2/3 DKO CD8 cells as the reference level
#Seperate up the sample_title column 
#Confirm that column names of count table and row names of the sample table are in the same order
sample_table <- sample_table %>% 
    separate(sample_title, c("egr_locus", "genotype", "cell_type", "replicates"), sep = "_") %>%
    unite(col = "condition", egr_locus, genotype, cell_type, sep = "_") %>% 
    dplyr::select(-c(species, library_layout)) %>% 
    mutate(condition=factor(condition, levels=c("Egr2/3_DKO_CD8", "Egr2/3_DKO_CD4", "Egr2_Kin_CD4", "Egr2_Kin_CD8")))
# View(sample_table)

levels(sample_table$condition)

#Make DESeq Data Set (dds) from matrix
dds <- DESeqDataSetFromMatrix(count_table, 
                              sample_table, 
                              ~condition)

colData(dds)
rowRanges(dds)

design(dds)
counts(dds) #equivalent to the command: assays(dds)$counts


#Access the design formula, counts matrix and sample informationfrom dds
#Calculate the size factors for each sample – estimateSizeFactors()
#Access the size factors from the dds object
#Generate a bar plot of the size factors for each sample, coloured by condition/group


dds <- estimateSizeFactors(dds)
sizeFactors(dds)

size_factors <- data.frame(sample = names(sizeFactors(dds)),
                           size_factor = sizeFactors(dds),
                           sample_group = colData(dds)$condition)


ggplot(size_factors, aes(y=size_factor, x=sample, fill=sample_group)) +
        geom_col() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Obtain dispersion estimates for each gene
dds <- estimateDispersions(dds)
dispersions(dds)
plotDispEsts(dds) 
#data is good quality, disp is affected by mean count, mean line central and final data forced towards fit line

#Perform the Wald test – nbinomWaldTest()
dds <- nbinomWaldTest(dds)


#Run all 3 functions (estimateSizeFactors(), estimateDispersions(), nbinomWaldTest())
dds <- DESeqDataSetFromMatrix(count_table, 
                              sample_table, 
                              ~condition)

dds_full <- DESeq(dds)

levels(sample_table$condition)
#Access the coefficients of the NB GLM - each coef is the baseline log expression compared against our comparitor (double KO CD8), then the next 3 columns is the change in the other conditions

head(coef(dds_full), 3)

#Access the results table for the comparison between CD8+ and CD4+ T cells from Egr2/3 DKO mice
res <- results(dds_full)
# View(res)


resultsNames(dds_full)
res <- results(dds_full, name = "condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8") #we came back and modified this after the MA plot


res_df <- as.data.frame(res)
# View(res_df)

# Plot a histogram of the raw and BH-adjusted p-values
# do they look as expected?

install.packages("patchwork")
install.packages("cowplot")

library(patchwork)
library(cowplot)


pval_plot <- ggplot(res_df, aes(x=pvalue))+
    geom_histogram()

padj_plot <- ggplot(res_df, aes(x=padj))+
    geom_histogram()


plot_grid(pval_plot, padj_plot, align = "h")


pval_plot + padj_plot

#cowplot - plot side by side

#reduced the no of pvals

#padj_plot / pval_plot #patchwork - plots in 1 col.


#Generate an MA plot of the log2 FC values for all genes
plotMA(res) #comparing CD4 to CD8 double KO

#Shrink the log2 FC values using the normal, apeglm and ashr methods

resNorm <- lfcShrink(dds_full, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = 'normal')

resAsh <- lfcShrink(dds_full, contrast = c('condition', 'Egr2/3_DKO_CD4', 'Egr2/3_DKO_CD8'), type = 'ashr')

resApeg <- lfcShrink(dds_full, coef = "condition_Egr2.3_DKO_CD4_vs_Egr2.3_DKO_CD8" , type = 'apeglm')


#compare methods with MA plots
plotMA(resAsh)
resAshplot <- recordPlot()

plotMA(resNorm)
resNormplot <- recordPlot()

plotMA(resApeg)
resApegplot <- recordPlot()

plot_grid(resAshplot, resNormplot, resApegplot, nrow = 3)
#note - used recordPlot to store plots produced in base R (MA plot)
#recommend using apglm as this is the default and tends to be appropriate for most data sets


#Generate a results table (one shrinkage method) containing mgi symbols
# Use the EnsDb.Mmusculus.v79 package
# Remove all genes with a padj of NA
# How many Ensembl IDs are not assigned an mgi symbol and how many mgi symbols are duplicated?

resApeg_table <- as.data.frame(resApeg)
# View(resApeg_table)

edb <- EnsDb.Mmusculus.v79
columns(edb)
gene_id <- select(edb, keys = rownames(resApeg_table), columns = c("GENEID", "SYMBOL"), keytype = "GENEID")
# View(gene_id)

gene_id$GENEID[duplicated(gene_id$GENEID)]
duplicated_genes <-  dplyr::filter(gene_id, SYMBOL %in% gene_id$SYMBOL[duplicated(gene_id$SYMBOL)])
# View(duplicated_genes)

resApeg_table <- rownames_to_column(resApeg_table, "GENEID") %>% 
    left_join(gene_id)
# View(resApeg_table)

sum(is.na(resApeg_table$SYMBOL))

resApeg_padj <- dplyr::filter(resApeg_table, !is.na(padj))
sum(is.na(resApeg_padj$padj))

#Write the results table to a CSV file
write.csv(resApeg_padj, file = "results/resApeg_padj.csv", quote = FALSE, row.names = FALSE)

#Filter the results table for padj < 0.05, and write to a CSV file
filtered_table <- dplyr::filter(resApeg_padj, padj < 0.05) %>% 
                    dplyr::filter(abs(log2FoldChange) > 1 )
# View(filtered_table)
dim(filtered_table)
write.csv(resApeg_padj, file = "results/filtered_table.csv", quote = FALSE, row.names = FALSE)

#Generate VST and rlog transformed counts:
#Plot the relationship between the mean expression and the sd of all genes – fit a trend line

vst <- vst(dds_full, blind = FALSE)
rlog <- rlog(dds_full, blind = FALSE)

assay(vst)

#Using both sets of transformed counts:
#Generate a PCA plot either using all genes, or top 500 most variable genes
#Generate a heatmap of the top 20 (by shrunken FC) differentially-expressed genes 
#label samples by condition and genes by mgi symbol
#Generate a heatmap of sample-sample correlations


vsn::meanSdPlot(assay(vst))

vsn::meanSdPlot(assay(rlog))
dim(vst)
plotPCA(vst, ntop = nrow(vst))
plotPCA(vst, ntop = 500)

# View(assay(vst))

top_20 <- filtered_table[order(-abs(filtered_table$log2FoldChange)), ] %>% 
    slice_head(n=20) %>% 
    pull(GENEID)
# View(top_20)
top_20

vst_top20 <- as.data.frame(assay(vst)) %>%
    dplyr::filter(row.names(.) %in% top_20)

pheatmap(vst_top20, scale = "row")


# volcano plot
# y = -log10(pvalue), x = fold change
df_for_volcanoplot <- 
    resApeg_padj %>% 
    mutate(y=-log10(padj)) %>% 
    mutate(highlight=factor(padj<0.05 & abs(log2FoldChange)>1))

color <- c('black', 'red')
color <- color[df_for_volcanoplot$highlight]

df_for_volcanoplot %>% 
    ggplot(aes(x=log2FoldChange, y=y)) +
    geom_point(color=color)

# input some random genes of interest (help, i have no idea??)
genes_of_interest <- c('Atp8b4', 'Ccdc7', 'Fscn1', 'Kcnj8', 'Rasd2', 'Sox8')

ggplot(df_for_volcanoplot, aes(x=log2FoldChange, y=y, label=SYMBOL)) +
    geom_point(color=color) +
    geom_label_repel(data = df_for_volcanoplot[df_for_volcanoplot$SYMBOL %in% genes_of_interest, ])




