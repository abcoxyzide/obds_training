lapply(c('tidyverse', 'cowplot', 'gridExtra','patchwork', 'scales'), require, character.only = TRUE)

#---
# Exercise
#---

## 0
# load data and add column
cds <- read.table('data/coding_gene_region.bed', col.names = c('chr', 'start','end','id', 'zero','strand'))
cds['width'] <- cds$end - cds$start

## 1
# Plot a histogram of the lengths using ggplot2:
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 10000)
# Add a plot title
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 10000) +
    labs(title = 'histogram')

# Change the x and y axis titles and sizes
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 10000) +
    labs(title = 'histogram',
         x = 'interval width',
         y = 'frequency') +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(face = 'italic'),
          plot.title = element_text(hjust = 0.5, size = 28))

# Change the size and angle of the x tick labels
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 10000) +
    labs(title = 'histogram',
         x = 'interval width',
         y = 'frequency') +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(face = 'italic'),
          plot.title = element_text(hjust = 0.5, size = 28)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(angle = 45, vjust = 0.5))

# Change the x axis upper limit to 500,000
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 10000) +
    labs(title = 'histogram',
         x = 'interval width',
         y = 'frequency') +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(face = 'italic'),
          plot.title = element_text(hjust = 0.5, size = 28)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(angle = 45, vjust = 0.5)) +
    # xlim(0, 500000) +
    scale_x_continuous(labels = comma, limits = c(0, 500000))

# Change the number of bins or the bin width
ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 2000) +
    labs(title = 'histogram',
         x = 'interval width',
         y = 'frequency') +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(face = 'italic'),
          plot.title = element_text(hjust = 0.5, size = 28)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(angle = 45, vjust = 0.5)) +
    # xlim(0, 500000) +
    scale_x_continuous(labels = comma, limits = c(0, 500000))

# Change the fill and border colour of the bars
p1 <- ggplot(cds, aes(x=width)) +
    geom_histogram(binwidth = 2000, color = 'red', fill = 'blue') +
    labs(title = 'histogram',
         x = 'interval width',
         y = 'frequency') +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(face = 'italic'),
          plot.title = element_text(hjust = 0.5, size = 28)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_text(angle = 45, vjust = 0.5)) +
    # xlim(0, 500000) +
    scale_x_continuous(labels = comma, limits = c(0, 500000))

## 3. 
# Using the diamonds dataset, plot a scatter plot of diamond length by price, coloured by the diamond colour and sized according to diamond width:
data(diamonds)

ggplot(data = diamonds, aes(x = x, y = price)) +
    geom_point(aes(size = y, colour = color))

# Use one of the ggplot built-in themes to alter the plot appearance
ggplot(data = diamonds, aes(x = x, y = price)) +
    geom_point(aes(size = y, colour = color)) +
    theme_minimal()

# Change the x axis upper limit to 12 and the intervals to 1.5
ggplot(data = diamonds, aes(x = x, y = price)) +
    geom_point(aes(size = y, colour = color)) +
    theme_minimal() +
    scale_x_continuous(breaks=seq(0,12,1.5))


# Add x and y axis titles and change their size
p2 <- ggplot(data = diamonds, aes(x = x, y = price)) +
    geom_point(aes(size = y, colour = color)) +
    theme_minimal() +
    scale_x_continuous(breaks=seq(0,12,1.5)) +
    labs(x = 'length',
         y = 'price') +
    theme(axis.title = element_text(size = 14))

# Plot the two plots that you have just made side-by-side using a ggplot2 function
#.. using patchwork library
p1 + p2
#.. using gridExtra
grid.arrange(p1, p2, ncol = 2, nrow = 1) # by default fills by row?
#.. using cowplot
plot_grid(p1, p2)
