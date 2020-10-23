library(Seurat)
library(tidyverse)
library(patchwork)
library(clustree)

pbmc.data <- Read10X(data.dir = "data/filtered_feature_bc_matrix/")

pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "pbmc3k")
pbmc[["ADT"]] <- CreateAssayObject(counts = pbmc.data$`Antibody Capture`)

Assays(pbmc)
DefaultAssay(pbmc) #RNA is the default already, don't need to change it yet

# nice command to view feature level metadata
View(pbmc[[]])

# continue
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # adding the percent.mt to the metadata slot as a column
pbmc[["percent.rib"]] <- PercentageFeatureSet(pbmc, pattern = "^RPS|^RPL")
# percent ribosomal gene sometimes used for filtering, eg if >50%; but need to be cautious

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA")) +
    geom_hline(yintercept=1250)
VlnPlot(pbmc, features = c("nCount_RNA"))
VlnPlot(pbmc, features = c("percent.mt")) +
    geom_hline(yintercept=14)
VlnPlot(pbmc, features = c("percent.rib"))




# visualize nCount RNA (ie UMI) and nFeature, colored by percent.mt (for QC)
pbmc[[]] %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    geom_point()
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# Suggestion: if just looking at this plot, can cut off at nFeature_RNA = 6000
# But this is an iterative process

# Common cutoffs:
# # 1000 for nFeature_RNA
# # 12.5 for percent.mt
# # Leave nCount_RNA for now
# # Leave percent.rib for now (but some people say that if it's above 50% in PBMCs then they should be filtered out, but be careful)

# Perform SCT
pbmc <- SCTransform(pbmc,
                    assay = 'RNA',
                    seed.use = 1448145,
                    verbose = TRUE)
# Remember that SCTransform performs these three functions:
# # NormalizeData
# # FindVariableFeatures
# # ScaleData

# The function creates a new assay called SCT, and changes the default to that
# common parameters to play around with:
# vars.to.regress ------ can regress out mt genes; BUT don't use it to regress batch effect, which has a huge systemic effect, which is better dealt with using methods like integration
# return.only.var.genes ------ reduces size of object and increase speed; but sometimes non-highly variable genes can be differentially expressed, so you might wanna tweak a bit once you have your clusters(?)

# Alternatively, can do the usual way of VST
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# see what highly variable genes are picked out (default number = 3000)
VariableFeatures(pbmc) %>% View

# scaling is done in sct
# run pca
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


# Visualize PC plot
DimPlot(pbmc, reduction = "pca")
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = percent.mt > 30)) # color by percent.mt

pbmc[['pca']]@cell.embeddings[, 1:2] %>%
    merge(pbmc[[]]['percent.mt'], by='row.names', all=T) %>% 
    ggplot(aes(x=PC_1, y=PC_2, color=percent.mt)) +
    geom_point()



# inspect PC loadings
print(pbmc[['pca']], dim=1:5, nfeatures=5)
# visualize
VizDimLoadings(pbmc, dims = 1, reduction = 'pca', balanced = F) # balanced defaulted to FALSE
DimHeatmap(pbmc, dims = 1, cells = 2000, balanced = T) # balanced = TRUE will give you equal number of over and underexpressed genes

# jackstraw doesn't work with SC transformed data
# go for elbow plot
ElbowPlot(pbmc, ndims = 50)
# in the vignette, it is said that we can be generous with SCT and pick a few more PC
# say 25


# Perform clustering
pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20) # num pc used is from above
# parameters to play with:
# # annoy.metric ------ default is euclidean, but cosine is better for some data
# # k.param ------ should play around with this iteratively, depending on size of clusters

# find cluster
pbmc <- FindClusters(pbmc, resolution = seq(0.2, 2, 0.2))
# depend on how many cells we have
# do this a lot of times, and later we'll check with clustree to find the optimal resolution (ie before it over-clusters)

# View(pbmc[[]])
# Switch to a resolution that you want to be in
Idents(pbmc) <- 'SCT_snn_res_0.8'


pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = percent.mt < 12.5))
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nCount_RNA > 20000))
DimPlot(pbmc, cells.highlight = WhichCells(pbmc, expression = nCount_RNA < 1000))
# you can see that all the dead cells / doublets cluster tgt, so very clean dataset

ggplot(pbmc[[]], aes(x=nFeature_RNA)) +
    geom_histogram(binwidth = 50, fill='grey', col='black')

# ---- Filtering ----
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & nFeature_RNA < 5500 & percent.mt < 12.5)
# then repeat the above steps
pbmc <- SCTransform(pbmc, assay = 'RNA', seed.use = 1448145, verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:25, k.param = 20)
pbmc <- FindClusters(pbmc, resolution = seq(0.2, 2, 0.2))

# run to detect any over clustering
clustree(pbmc)
# seems to be good until 1.4-1.6
# but we're not interested in b cell subclusters, so let's pick 0.8
Idents(pbmc) <- "SCT_snn_res.0.8"

pbmc <- RunUMAP(pbmc, dims = 1:25)
DimPlot(pbmc, reduction = "umap")
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# ---- ADT ----
DefaultAssay(pbmc) <- 'ADT'


# see what antibody we have
# (if we were in RNA, then we'd see all the genes)
rownames(pbmc) %>% View

# recommended normalization method for ADT data; may change in future
pbmc <- NormalizeData(pbmc, assay = 'ADT', normalization.method = 'CLR')
pbmc <- ScaleData(pbmc, assay = 'ADT')

## visualize antibody on clusters based on RNA
DefaultAssay(pbmc) <- 'RNA' # weird error, can't find gene in ADT, this is a workaround
# color the UMAP
FeaturePlot(pbmc, features = c("CD3-TotalSeqB", "CD19-TotalSeqB", "CD4-TotalSeqB", "CD8a-TotalSeqB", "CD14-TotalSeqB",
                               "CD3E", "CD19", "CD4", "CD8A", "CD14"), 
            min.cutoff = "q05", max.cutoff = "q95", ncol = 5)
# ridgeplot
RidgePlot(pbmc, features = c("CD19-TotalSeqB", "CD45RA-TotalSeqB", "CD20-TotalSeqB", "PD-1-TotalSeqB"), ncol = 2)
# revert to ADT
DefaultAssay(pbmc) <- 'ADT'



# get rid of the IgG controls
features_adt <- grep(pattern = 'control', rownames(pbmc), invert = T, value = T)

# Run PCA
pbmc <- RunPCA(pbmc,
               features = features_adt,
               npcs = 14,
               reduction.name = 'pca_adt',
               reduction.key = 'adtPC_'
               )
Reductions(pbmc)
DimPlot(pbmc, reduction = 'pca_adt')
ElbowPlot(pbmc, ndims = 14, reduction = 'pca_adt')
# 10 seems to be reasonable


# run UMAP, while still colored by clustering based on RNA
pbmc <- RunUMAP(pbmc, dims = 1:10,
                reduction = 'pca_adt',
                reduction.name = 'umap_adt',
                reduction.key = 'adtUMAP_')
DimPlot(pbmc, reduction = 'umap_adt')




# now run clustering based on adt
# before that, let's stash the active ident
pbmc[['rnaCluster']] <- Idents(pbmc)

pbmc <- FindNeighbors(pbmc, k.param = 20, reduction = 'pca_adt', dims = 1:10)
# alternatively, we can directly use the proteins if the dimension isn't high
# pbmc <- FindNeighbors(pbmc, k.param = 20, features = features_adt, dims = NULL)
pbmc <- FindClusters(pbmc, resolution = seq(0.1, 0.8, 0.1), graph.name = "ADT_snn")

clustree(pbmc)
# 0.4 seems ideal
Idents(pbmc)<-"ADT_snn_res.0.4"

# run and plot UMAP again with res 0.4
pbmc <- RunUMAP(pbmc, dims = 1:10,
                reduction = 'pca_adt',
                reduction.name = 'umap_adt',
                reduction.key = 'adtUMAP_')
DimPlot(pbmc, reduction = 'umap_adt')

# stash the cluster ID
pbmc[['adtCluster']] <- Idents(pbmc)

clustering.table <- table(pbmc[[]][['adtCluster']], pbmc[[]]$SCT_snn_res.0.8)

