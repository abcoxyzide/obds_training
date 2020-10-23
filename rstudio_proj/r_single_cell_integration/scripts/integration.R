library(scran)
library(scater)
library(DropletUtils)
library(Seurat)
library(batchelor)
library(tidyverse)

sce <- read10xCounts(c(v2 = 'data/pbmc_1k_v2/',
                       v3 = 'data/pbmc_1k_v3/'))

# using pre-filtered data, already got rid of empty droplets
# quick standard filtering
is.mito <- grepl("^MT-", rowData(sce)$Symbol)
table(is.mito)
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
# function:     - uses isOutlier function to classify the cells; caution that the 'MAD' (median absolute deviation) statistic is used, which can be different depending on the direction you're looking at (higher / lower)
table(filtered)
sce <- sce[, !filtered$discard]

### Pre-integration visualization ----
# go straight to visualization to see if integration is needed between the two datasets with different chemistry
sce <- logNormCounts(sce)
allf <- modelGeneVar(sce)
hvg <- getTopHVGs(allf, var.field = 'bio', n = 1000) # save this variable, so we can use it down there for fair comparison pre- and post-integration
rowData(sce)['HVG'] <- rownames(sce) %in% hvg

set.seed(123)
sce <- runPCA(sce, ntop = Inf, subset_row = rowData(sce)[['HVG']])

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")

plotReducedDim(sce, 'PCA', colour_by = 'Sample')
sce <- runUMAP(sce, dimred = 'PCA', n_dimred = 10)
plotReducedDim(sce, dimred = 'UMAP', colour_by = 'Sample')

# Comment:
# there's a systemic shift to the left for v3 chemistry
# and it looks horrible on UMAP


### Integration ----

# first split the data into two SCE objects
sce_v2 <- sce[, sce$Sample=='v2']
sce_v3 <- sce[, sce$Sample=='v3']
