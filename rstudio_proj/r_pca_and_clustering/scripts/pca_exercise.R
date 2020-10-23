sapply(c('tidyverse', 'cowplot', 'umap'), require, character.only=T)

df_logcounts <- read.csv('data/logcounts.csv', row.names = 1)
df_cellmetadata <- read.csv('data/cell_metadata.csv', row.names=1)
df_genemetadata <- read.csv('data/gene_metadata.csv', row.names=1)

df_cellmetadata$Lane <- factor(df_cellmetadata$Lane)

#-----
# do pca on unscaled count table data
pc_logcounts <- df_logcounts %>%
    as.matrix %>% t %>%
    prcomp

# create dataframe
df_pc_var <- tibble(
    PC = paste0('PC', seq_along(pc_logcounts$sdev)),
    explained_var = pc_logcounts$sdev^2/sum(pc_logcounts$sdev^2)*100,
    cumulative_var = cumsum(explained_var)
) %>%
    mutate(PC = factor(PC, PC)) # necessary to set the levels as is, so that R doesn't arrange them by alphabetical order (ie PC1, PC10, PC100, etc)

# standard scree plot
df_pc_var %>%
    ggplot(aes(x=PC, y=explained_var)) +
    geom_point()

# variant of scree plot, based on cumulative variance
df_pc_var %>%
    ggplot(aes(x=1:length(explained_var), y=cumulative_var)) +
    geom_point()

#-----
# build dataframe of PCA outputs to be tidy + has metadata
df_pctransformed <- pc_logcounts$x %>% as.data.frame %>%
    rownames_to_column('cellid') %>%
    left_join(rownames_to_column(df_cellmetadata, 'cellid'))


# visualize
plist <- lapply(
    colnames(df_cellmetadata),
    function(col) {
        df_pctransformed %>% as_tibble %>%
            mutate(color=.[[col]]) %>%
            ggplot(aes(x=PC1, y=PC2, color=color)) +
            geom_point() +
            theme(legend.position="none")
    }
)
plot_grid(plotlist = plist, ncol=3)

#-----
# the "rotation" attr is the pc loadings (ie significance / eigenvectors)




#-----
# kmeans

km <- lapply(2:10, function(i) kmeans(t(df_logcounts), centers = i)) # centers is just the number of cluster u think there are



df_pctransformed$cluster <- as.factor(km[[3]]$cluster[df_pctransformed$cellid])
ggplot(df_pctransformed, aes(x=PC1, y=PC2, color = cluster)) +
    geom_point()

df_km_out <- sapply(
    2:10,
    function(i) {
        km_object <- kmeans(t(df_logcounts), centers = i)
        withinss <- sum(km_object$withinss)
        b_t <- km_object$betweenss / km_object$totss
        return(c(i, withinss, b_t))
    }
) %>% t %>%
    as.data.frame %>%
    set_names(c('k', 'withinss', 'b_t'))

#
df_km_out %>%
    pivot_longer(-k) %>%
    ggplot(aes(x=k, y=value, color=name))+
    geom_line()

###

umap_object <- umap(pc_logcounts$x[, 1:12])
head(umap_object$layout)
umap_coord <- as.data.frame(umap_object$layout) %>%
    set_names(c('umap1', 'umap2')) %>% rownames_to_column('cellid')

umap_coord %>%
    left_join(rownames_to_column(df_cellmetadata, 'cellid')) %>%
    ggplot(aes(x=umap1, y=umap2, color=Time)) +
    geom_point() +
    facet_wrap(~Infection)
