###22/10 exercise scRNA using bioconductor###

library(tidyverse)
library(scater)
library(scran)
library(DropletUtils)
library(DelayedMatrixStats)
library(DelayedArray)
library(Matrix)
library(sparseMatrixStats)


#importing 'filtered' matrix using DropUtils
sc_data <- read10xCounts(c(filtered = "bioc/filtered_feature_bc_matrix"), col.names = TRUE) #True allows you to see colnames rather than null.
sc_data #33538 x 1221
slotNames(sc_data) #indicates what slots you have

colnames(colData(sc_data)) #sample & barcode
colData(sc_data) #table format of above info #c(filtered = "file path") from line 10 will have less info under sample col (just filtered)

colnames(rowData(sc_data)) #ID, symbol and type (gene or protein)
rowData(sc_data) #gives some table content of above.

assayNames(sc_data) #counts

metadata(sc_data) #tells you nickname and file path

assay(sc_data)

barcode <- barcodeRanks(assay(sc_data)) #rank, total and fitted (all numeric) > this is used to build the plot in p14 slide
#if dont specificy which assay it will pick up the 1st so to be specific can use following codes:
    #assay(sc_data, "counts") or assays(sc_data)$counts
class(barcode) #Dframe
ggplot(data.frame(barcode), aes(x=rank, y=total, col=fitted)) +
    geom_point() +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(trans = "log10") #or use scale_x_log10() as shorthand.

#compute/visualize QC metrix using scater to remove any cells?
#per cell QC
per_cell <- perCellQCMetrics(sc_data) #Compute per-cell quality control metrics for a count matrix or a SingleCellExperiment.
#argument ; %_top (% UMI accounted for by the top n genes, want to capture a diverse transcriptome), exp vals (specify name of count matrix, "counts" by default)
per_cell #DataFrame, check under help (value)
#eg. first cell has 8512 UMI/ncounts, distributed to 2663 genes.

ggplot(data.frame(per_cell), aes(x = percent_top_50)) +
    geom_histogram(color = "black", fill = "white", bins = 50) #counting per cell. each cell what fraction of UMI in that cell would account for. In this case, >60% of library poor quality where few transcripts were captured (but may be another population of cells)
#%_top = % of counts assigned to %_toppage of most highly expression genes.

#scatter plot
ggplot(data.frame(per_cell), aes(x = percent_top_50, y = total)) +
    geom_point() #low total count and high percent_top_50 so want to remove them
#correlation high val vs small library size
#could plot against mt to see if these in the bottom right are dying cells.

#trying out another strategy
sc_data <- addPerCellQC(sc_data) #compute QC metrics and add them to a SummarizedExperiment's metadata.
sc_data #overwritting
colData(sc_data)

ggplot(as.data.frame(colData(sc_data)), aes(x = percent_top_50)) +
    geom_histogram(color = "black", fill = "white", bins = 50)


#calculate per feature QC
per_feature <- perFeatureQCMetrics(sc_data)
per_feature #mean count/feature & detected (fraction of samples the genes have been detected)

sc_data <- addPerFeatureQC(sc_data) #'stuffing' to same obj
rowData(sc_data) #added mean & detected (n=5cols)

#we are NOT filtering yet! ---- collecting 'evidence' first (some cells may looked bad intially but arent)

#convert counts > normalised exp vals to remove cell-specific biases (scater/scran)

sc_data <- logNormCounts(sc_data) #calculating both size factor % normliased log counts.
#adjusting for library size and normalising it > work on log_CPM. usually counts/10,000 in scRNA.
sc_data
range(assay(sc_data, "logcounts")) #0.00000-12.68275

#
rowMeans(assay(sc_data, "logcounts")) #need to specify 'assay' otherwise will pick the 1st (count)
sparseMatrixStats::rowVars(assay(sc_data, "logcounts")) #est. variance by row while minimising memory

mean_var <- data.frame(
    mean = rowMeans(assay(sc_data, "logcounts")),
    variance = sparseMatrixStats::rowVars(assay(sc_data, "logcounts"))
)
#mean as 1st col;
#mean on x, var on y
ggplot(data.frame(mean_var), aes(x = mean, y = variance)) +
    geom_point() #post-normalisation

#plotting on raw data (against normalisation)
mean_var <- data.frame(
    mean = rowMeans(log(assay(sc_data, "counts") + 1)),
    variance = rowVars(log(as.matrix(assay(sc_data, "counts") + 1)))
    )

ggplot(data.frame(mean_var), aes(x = mean, y = variance)) +
    geom_point() #lower mean and var pre-normalization


#filtering-select features of downstream analyses e.g. HVGs using scran.

decomposed_log <- modelGeneVar(sc_data) #modelling var of log-exp for each gene based on mean-Var trend.
decomposed_log #technical/biolgoical component of the var.
#apply DR to compact data and further reduce noise using scater

plot(decomposed_log$mean, decomposed_log$total)
curve(metadata(decomposed_log)$trend(x), add=TRUE, col="dodgerblue")
#each dot = gene; line = trend of total var (generated where most genes are)
#HVGs e.g. the top 3 dots > driving heterogeneity but want to grab enough genes (above blue line)
#biological variance is the dots above the line (beyond technical); technical var = line.

#define HVGs (above blue line)
getTopHVGs(decomposed_log,
           var.field = "bio",
           n = 1000)

hvgs <- getTopHVGs(decomposed_log,
                   var.field = "bio",
                   n = 1000)

#Apply dimensionality reduction to compact the data and further reduce noise
set.seed(123)

sce <- runPCA(sce, name = "PCA",
              ncomponents = 50,
              ntop = Inf,     #use all the genes
              subset_row = hvgs,
              scale = FALSE
)

reducedDims(sce)
str(reducedDims(sce)$PCA)

plotReducedDim(sce,
               dimred = "PCA",
               colour_by = "ENSG00000148737") #colData(sce)  colour_by = "percent_top_50"

percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab="PC", ylab="Variance explained (%)")

sce <- runUMAP(sce,
               name = "UMAP",
               dimred = "PCA",
               n_dimred = 30) #choose the top 30 components

plotReducedDim(sce,
               dimred = "UMAP",
               colour_by = "percent_top_50")

hvgs

SNNGraph <- buildSNNGraph(sce, use.dimred = 'PCA')
colData(sce)[["cluster_louvain"]] <- factor(igraph::cluster_louvain(SNNGraph)$membership)

colData(sce)[["cluster_louvain"]]

plotReducedDim(sce,
               dimred = "UMAP",
               colour_by = "cluster_louvain")

#identify markers per cluster using scran
markers <- findMarkers(sce, groups = sce$cluster_louvain)
markers[[1]] #cluster compared with other clusters

library(iSEE)

iSEE(sce)
