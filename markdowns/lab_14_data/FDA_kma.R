# An Introduction to Functional Data Analysis
# Applied Statistics, 2021/2022 
# 

# 
# K-mean alignment
# 

setwd("~/Corsi/Statistica Applicata/Applied Statistics MATE 21-22/Lab 14 - 31052022")

# Based on Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. 
# "K-mean alignment for curve clustering". 
# Computational Statistics and Data Analysis, 54, 1219-1233.


# Using the R package fdakma
library(fdakma)
help(kma)

help(kma.data)

data(kma.data)
names(kma.data)

x <- kma.data$x   # abscissas
y0 <- kma.data$y0 # evaluations of original functions
y1 <- kma.data$y1 # evaluations of original functions' first derivatives

# Plot of original functions
x11()
matplot(t(x),t(y0), type='l', xlab='x', ylab='orig.func')
title ('Original functions')

# There seems to be three clusters of functions obtained by warping the 
# abscissas.

# Plot of original function first derivatives
x11()
matplot(t(x),t(y1), type='l', xlab='x', ylab='orig.deriv')
title ('Original function first derivatives')


# Without alignment, let's try with 3 clusters

set.seed(4)
fdakma_example_noalign_0der <- kma(
  x=x, y0=y0, n.clust = 3, 
  warping.method = 'NOalignment', 
  similarity.method = 'd0.pearson',   # similarity computed as the cosine
                                      # between the original curves 
                                      # (correlation)
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example_noalign_0der)
fdakma_example_noalign_0der$labels

# We might want to use a similarity measure which considers first 
# derivatives

set.seed(5)
fdakma_example_noalign_1der <- kma(
  x=x, y0=y0, y1=y1, n.clust = 3, 
  warping.method = 'NOalignment', 
  similarity.method = 'd1.pearson',   # similarity computed as the cosine
                                      # between the first derivatives 
                                      # (correlation)
  center.method = 'k-means'
  #,seeds = c(1,11,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example_noalign_1der)
fdakma_example_noalign_1der$labels

# Label switching aside, it is the same clustering

table(fdakma_example_noalign_0der$labels,
      fdakma_example_noalign_1der$labels, dnn=c("0der", "1der"))


# Let's allow for alignment. 
# Result of kma with 2 clusters, 
# allowing affine transformation for the abscissas
# and considering 'd1.pearson' as similarity.method.

set.seed(1)
fdakma_example <- kma(
  x=x, y0=y0, y1=y1, n.clust = 2, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',  # similarity computed as the cosine
                                     # between the first derivatives 
                                     # (correlation)
  center.method = 'k-means'
  #seeds = c(1,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)

# Much better: high within group similarity after alignment, 
# the boxplot is concentrated on 1.

# Labels assigned to each function
fdakma_example$labels

table(fdakma_example_noalign_0der$labels,
      fdakma_example$labels, 
      dnn=c("NOalign_0der_3groups", "Align_1der_2groups"))

# Total shifts and dilations applied to the original 
# abscissa to obtain the aligned abscissa
fdakma_example$shift
fdakma_example$dilation


# Warpings applied to the original data, colored according to  
# membership to the 3 clusters obtained via k-means, no alignment,
# d0.pearson similarity

plot(x, type="n", xlim=c(min(x),max(x)), ylim=c(min(x),max(x)+2), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:30)(
  abline(a=fdakma_example$shift[i],b=fdakma_example$dilation[i], 
         col=fdakma_example_noalign_0der$labels[i])
)  


# How to choose the number of clusters and the warping method
help(kma.compare)

kma.compare_example <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=TRUE)

kma.compare_example_2 <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c("shift", "dilation", "affine"), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=TRUE)

kma.compare_example_3 <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c("NOalignment", "shift", "dilation", "affine"), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=TRUE)
