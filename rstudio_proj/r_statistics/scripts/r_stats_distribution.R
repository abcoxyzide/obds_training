sapply(c('tidyverse'), require, character.only=T)


# Set seed as today's date! :D
set.seed(20201014)

# make table with various distribution
dst <- tibble(rnorm=rnorm(1000))

#hist(dst$rnorm, breaks = 25)

# practice ploting using ggplot
dst %>%
    ggplot(aes(x=rnorm)) +
    geom_histogram() +
    geom_vline(xintercept = 10)

#-----------------------------------------#
# For the standard normal distribution (mean = 0, sd = 1)
# 1. Plot the cumulative distribution function in the range [-5, 5]
# 2. Plot the inverse cumulative distribution function for quantiles in 0.01 increment.
# 3. Plot the density function in the range [-5, 5]
# 4. What is the probability of observing a value greater than 2?
# 5. What is the probability of observing a value between -2 and 2?
# 6. What is the probability of observing a value more extreme than -2 or 2?
#-----------------------------------------#

## /Q1/
# set up the x coordinate
x_index = seq(-5, 5, length.out = 1000)

# generate the required function of the distribution, and input it as column in a tidyverse table (ie tibble)
dst <- tibble(
    x_coordinate = x_index,
    cdf_normal_distribution = pnorm(q = x_index, mean = 0, sd = 1),
)

# plot using ggplot
dst %>%
    ggplot(aes(x=x_coordinate, y=cdf_normal_distribution)) +
    geom_point()

## /Q2/
# set up the x coordinate for probability
x_prob = seq(0, 1, 0.01)# length.out = 1000)

dst <- tibble(
    x_prob = x_prob,
    inverse_cdf_normal = qnorm(p = x_prob, mean=0, sd=1)
)

dst %>%
    ggplot(aes(x=x_prob, y=inverse_cdf_normal)) +
    geom_point()

## /Q3/
# set up the x coordinate
x_index = seq(-5, 5, length.out = 1000)
# generate the required function of the distribution, and input it as column in a tidyverse table (ie tibble)
dst <- tibble(
    x_coordinate = x_index,
    density_normal_distribution = dnorm(x = x_index, mean = 0, sd = 1)
)
# plot using ggplot
dst %>%
    ggplot(aes(x=x_coordinate, y=density_normal_distribution)) +
    geom_point()

## /Q4/
1 - pnorm(2)

## /Q5/
pnorm(2) - pnorm(-2)

## /Q6/
1 - (pnorm(2) - pnorm(-2))

#-----------------------------------------#
# Testing out the ecdf and knots function #
#-----------------------------------------#
ecdf_object = ecdf( iris$Sepal.Length )
knots_object = knots(ecdf_object)
plot(ecdf_object)
plot(knots_object)

#---------------#
# more exercise #
#---------------#
# Visualise the distribution of Sepal.Length , stratified by species #
iris %>%
    ggplot(aes(Sepal.Length)) +
    geom_histogram(color='black') +
    facet_wrap(~Species, ncol=1)

# are they normally distributed?
shapiro.test(iris$Sepal.Length)
# no not as a whole
# per species?
iris %>%
    group_by(Species) %>%
    summarise(shapiro=shapiro.test(Sepal.Length)$p.value)
# yes they are within each species!

# significant variation of Sepal.Length between the various species?
# do one way anova (>2 samples, normal dist, unpaired)
aov(formula = Sepal.Length ~ Species, data = iris) %>%
    summary

#---------------#
# more exercise #
#---------------#
data(ChickWeight)
lapply(ChickWeight, summary)

shapiro.test(ChickWeight$weight)

# univariate linear model on (weight ~ Diet), and (weight ~ Time)
# look up ?formula for more info; other operators are : ^ *
# the model (weight ~ Diet) is using Diet 1 as reference, essentially the first guy in the levels(data); you can change that with relevel(data, 'i_want_this_first')


formula_1 <- formula(weight~Diet+Time)
lm_chick <- lm(formula = formula_1, data = ChickWeight)
summary(lm_chick)

# multivariate linear model
formula_2 <- formula(weight~Diet*Time)
lm_chick <- lm(formula = formula_2, data = ChickWeight)
summary(lm_chick)


ChickWeight %>%
    ggplot(aes(x = Time, y = weight, color = Diet)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    geom_abline(slope = 6.8418, intercept = 30.9310) + # the line for Diet1:Time (ie Time)
    geom_abline(slope = 6.8418 + 4.5811, intercept = 30.9310 - 12.6807) # the line for Diet3:Time


#---------------#
# more exercise #
#---------------#
mtx_count <- read.csv('data/logcounts.csv', row.names = 1)
mtx_meta <- read.csv('data/cell_metadata.csv')

tdy_data <- mtx_count %>%
    t %>% as.data.frame %>%
    rownames_to_column('cellid') %>%
    left_join(mtx_meta, by=c('cellid'='Sample')) %>% tibble

# perform t test for every gene
p_values <- vapply(select(tdy_data,-c(Infection, cellid)), function(x) wilcox.test(x ~ tdy_data$Infection)$p.value, 0)
hist(p_values, breaks=100)

# then adj the p values for multiple testing

test <- lapply(p.adjust.methods, function(padj_method) {

padj_values <- p.adjust(p_values, method=padj_method)
hist(padj_values, breaks=100)

#=====================#
'
Turns out that the method of p value correction can vastly affect the results
Might be due to multiple factors, eg using t test instead of wilcoxin, or that
not the entire genome is given
'
#=====================#

padj_tbl <- data.frame(padj=padj_values) %>%
    rownames_to_column('geneid') %>%
    mutate(sig=(padj<0.05))

#-----------------
# over representation analysis
#-----------------
go_info <- read.csv('data/go_info.csv')
human_go <- read.csv('data/human_go_bp.csv')

list_go <- split(human_go$ensembl_gene_id, human_go$go_id)

list_contingency_mtx <-
with(padj_tbl,
lapply(list_go, function(go_term){
    sigf <- table(factor(geneid[sig] %in% go_term, levels = c(T, F))) # by specifying this vector as factor and setting the corresponding values, we can force the table to have the corresponding outcomes; without which is gonna break the code, cus it'll return only 1 column sometimes
    isgf <- table(factor(geneid[!sig] %in% go_term, levels=c(T, F)))
    rbind(sigf,isgf)
}))

p_fisher <-
    sapply(list_contingency_mtx, function(m) {
    fisher.test(m, alternative = 'greater')$p.value
})

df_ora <-
    data.frame(pfisher=p_fisher) %>% rownames_to_column('go_id') %>%
    mutate(padj = p.adjust(pfisher, method=padj_method)) %>%
    left_join(go_info, by=c('go_id'='GOID')) %>%
    arrange(padj)

return(head(df_ora, n=20))
})
