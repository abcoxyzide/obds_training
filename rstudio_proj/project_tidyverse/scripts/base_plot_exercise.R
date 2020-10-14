# renv::init()
# install.packages("BiocManager")
# BiocManager::install(c('tidyverse', 'cowplot', 'gridExtra','patchwork', 'biomaRt', 'pheatmap'))

lapply(c('tidyverse', 'cowplot', 'gridExtra','patchwork'), require, character.only = TRUE)

#---
# Exercise
#---

## 1
# Draw a bar plot showing the number of cars with 3, 4 or 5 forward gears. Hint: use table() function to get data into the correct format.
data(mtcars)
barplot(table(mtcars$gear))

# Change plot margins to 6, 6, 5, 5
p = par(no.readonly = T)
par(mar=c(6,6,5,5))
barplot(table(mtcars$gear))

# Add x and y axis labels
barplot(table(mtcars$gear),
        xlab = 'num of gears',
        ylab = 'num of cars')

# Add a main title spread across two lines
barplot(table(mtcars$gear),
        xlab = 'num of gears',
        ylab = 'num of cars',
        main = 'Barplot of number of cars\nwith number of gears')

# Change the colour and width of the bars
barplot(table(mtcars$gear),
        xlab = 'num of gears',
        ylab = 'num of cars',
        main = 'Barplot of number of cars\nwith number of gears',
        col = 'white', border = 'red',
        width = c(1,4,10), space = c(0, 1, 3)
        )

# Add a horizontal line at y = 6
abline(h=6)

## 2.
# Generate a scatter plot of mpg vs. hp coloured by gear values
par(p)
with(mtcars, plot(mpg, hp, col=gear))
# Change points to filled circles and increase their size
with(mtcars, plot(mpg, hp, col=gear,
                  pch=16, cex=0.5))
# Add x and y axis titles and change size
with(mtcars, plot(mpg, hp, col=gear,
                  pch=16, cex=0.5,
                  xlab = 'mpg of mtcars',
                  ylab = 'hp of mtcars',
                  cex.lab = 1.3))
# Add a legend
legend(legend = unique(mtcars$gear), col=unique(mtcars$gear), pch=16, x = 'right')
# Change colours to red, green and blue
with(mtcars, plot(mpg, hp, col=c('red', 'green', 'blue')[factor(gear)],
                  pch=16, cex=0.5,
                  xlab = 'mpg of mtcars',
                  ylab = 'hp of mtcars',
                  cex.lab = 1.3))
legend(legend = levels(factor(mtcars$gear)), col=c('red', 'green', 'blue'), pch=16, x = 'right')

## 3.
# Plot the two plots that you have just made side-by-side by changing the global graphical parameters
par(mfrow=c(1,2))
barplot(table(mtcars$gear),
        xlab = 'num of gears',
        ylab = 'num of cars',
        main = 'Barplot of number of cars\nwith number of gears',
        col = 'white', border = 'red',
        width = c(1,4,10), space = c(0, 1, 3)
)
abline(h=6)
with(mtcars, plot(mpg, hp, col=c('red', 'green', 'blue')[factor(gear)],
                  pch=16, cex=0.5,
                  xlab = 'mpg of mtcars',
                  ylab = 'hp of mtcars',
                  cex.lab = 1.3))
legend(legend = levels(factor(mtcars$gear)), col=c('red', 'green', 'blue'), pch=16, x = 'right')
