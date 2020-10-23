library(tidyverse)
library(magrittr)
library(scater)
library(scran)
library(DropletUtils)
library(PCAtools)
library(iSEE)


## Pre-processing -----
sce <- read10xCounts(c(filtered = 'data/filtered_feature_bc_matrix/')) # doing so will name your sample as 'filtered'
sce

# Basic operations to access what is in the data object
slotNames(sce) # equivalent to str() but only printing the names so readable
colData(sce)
rowData(sce)
assayNames(sce)
assay(sce)[1:6, 1:36] # note there is no col names; to have it, turn on col.name=T when reading in the data
# assay(sce, 'count') # to explicitly access the assay, in case there are more than one of them; same logic for reducedDim()
# retrieveCellInfo() is also potentially useful

metadata(sce)



barcode <- barcodeRanks(assay(sce))
barcode %>%
    as.data.frame %>%
    ggplot(aes(x=rank, y=total, col=fitted)) +
    geom_point() +
    scale_y_continuous(trans = 'log10') +
    scale_x_log10() # shorthand of the above line for x axis


# cell QC as a separate object
percellQC <- perCellQCMetrics(sce, exprs_values = "counts",
                              percent_top = seq(1, 101, 10)) %T>% print
# option:   subset        - can calculate QC metrics also on the subset you specified
# option:   percent_top   - the percent of UMI that the top N genes account for
# output:   sum           - number of UMI
# output:   detected      - number of genes that has count > detection_limit (default = 0)

# but standard practice is to add it inside the singlecellexperiment
sce <- addPerCellQC(sce, exprs_values = "counts",
                    percent_top = seq(20, 200, 20))
colData(sce)

colData(sce) %>% as.data.frame %>%
    pivot_longer(-c(Sample, Barcode, sum, detected, total)) %>%
    mutate(name = factor(name, levels = unique(name))) %>%
    ggplot(aes(x=value, y = total)) +
    geom_point() +
    facet_wrap(~name)


# there's also feature QC
perFeatureQCMetrics(sce, exprs_values = 'counts')
# output:   mean        - average count of that gene
# output:   detected    - the percentage of observations > detection_limit (default = 0)

sce <- addPerFeatureQC(sce, exprs_values = 'counts')
rowData(sce)

rowData(sce) %>% as.data.frame %>%
    ggplot(aes(x = mean)) +
    geom_histogram(bins = 50)

# leave them for now, make graphs first and filter later iteratively

## Normalization -----

sce <- logNormCounts(sce)
# // Function:
#   - make a new assay called logcounts
#   - calculate a quick sizefactor per cell

range(assay(sce, 'logcounts'))
range(assay(sce, 'counts'))


# make some plots to show that the normalization did
# before normalization, counts would be affected by library size
g <- assay(sce, 'counts') # raw counts
tail(sort(rowMeans(g))) # then cherry picked on one of the more expressed genes
g <- g['ENSG00000198938',]
libsize <- colData(sce)[['total']]
data.frame(
    gene = g,
    libsize = libsize
) %>%
    ggplot(aes(x=gene, y=libsize)) +
    geom_point()

g2 <- assay(sce, 'logcounts')['ENSG00000198938',]
g2 <- 2^g2 # reverse the log2 transform, so now the effect would be due to normalization
data.frame(
    gene_normalized = g2,
    libsize = libsize
) %>%
    ggplot(aes(x=gene_normalized, y=libsize)) +
    geom_point()
# sparse matrix is memory efficient; but needs special packages to handle if you don't want to convert it into a regular matrix (all the 0 will now take up space)


## Feature selection -----

# calculate the per-gene variance
allf <- modelGeneVar(sce) %T>% print

# visualizing, where the blue line is the overall variance
plot(allf$mean, allf$total, col=rgb(0,0,0, 0.5))
curve(metadata(allf)$trend(x), add=TRUE, col="dodgerblue")
# plot: - each dot = 1 gene
# plot: - blue line = where most of the genes are
# plot: - HVG (highly variable gene) eg the top 3, driving heterogeneity

hvg <- getTopHVGs(allf,
           var.field = 'bio', # the 'real' biological difference, as estimated by the deconvolution method in modelGeneVar
           n = 1000) # pick a number high enough to tell difference, but small enough to run smoothly & reduce multiple testing

# add to the gene metadata
rowData(sce)['HVG'] <- rownames(sce) %in% hvg

set.seed(123) # some stochastic processes are used in this package
sce <- runPCA(sce,
              ntop = Inf,
              subset_row = rowData(sce)[['HVG']],
              scale = F)
# method:       - runPCA returns SCE object; calculatePCA returns a PCA object
# ncomponenets: - return only 50 PCs by default
# ntop:         - using 500 most variable genes max by default
# scale:        - FALSE by default

# access the PCA object generated above
summary(reducedDim(sce, 'PCA'))

plotReducedDim(sce, dimred = 'PCA', ncomponents = 3)
plotReducedDim(sce, dimred = 'PCA', ncomponents = 2)

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
abline(v = chosen.elbow)


# now run a UMAP
sce <- runUMAP(sce,
        dimred = 'PCA',
        n_dimred = 6,
        name = 'UMAP')
plotReducedDim(sce, dimred = 'UMAP',
               colour_by = 'percent_top_100')
plotReducedDim(sce, dimred = 'UMAP',
               colour_by = 'ENSG00000198938')


# trying out the quickcluster function
# usually with more noise
colData(sce)$quickCluster <- quickCluster(sce, assay.type="logcounts")
plotReducedDim(sce,
               dimred = "UMAP",
               ncomponents = 2,
               colour_by="quickCluster")


# cluster by PCA
snn <- buildSNNGraph(sce, use.dimred = 'PCA')
colData(sce)[["cluster_louvain"]] <- factor(igraph::cluster_louvain(snn)$membership)
plotReducedDim(sce, dimred = 'UMAP',
               colour_by = 'cluster_louvain')

# cluster by UMAP? is that valid?
snn <- buildSNNGraph(sce, use.dimred = 'UMAP')
colData(sce)[["cluster_louvain_UMAP"]] <- factor(igraph::cluster_louvain(snn)$membership)
plotReducedDim(sce, dimred = 'UMAP',
               colour_by = 'cluster_louvain_UMAP')


markers <- findMarkers(sce, groups = sce$cluster_louvain)

markers[[1]] #cluster 1 compared with other clusters; need to do some calculation / merge the results to get the overall marker genes?

iSEE(sce)
















