library(fdakma)

################### K-MEAN ALIGNMENT #####################
#---------------------------------------------------------
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



#### WITHOUT ALIGNMENT
###------------------
# Try 3 clusters n.clust = 3

# 1) without first derivative
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


# 2) with first derivative
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


# Compare subdivision in groups
table(fdakma_example_noalign_0der$labels,
      fdakma_example_noalign_1der$labels, dnn=c("0der", "1der"))






#### WITH ALIGNMENT
###------------------
# Try with 2 clusters

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

# Warpings applied to the original data, colored according to  
# membership to the 3 clusters obtained via k-means, no alignment,
# d0.pearson similarity
plot(x, type="n", xlim=c(min(x),max(x)), ylim=c(min(x),max(x)+2), xlab="abscissa", ylab="warping")
title("Alignment affinities")
for(i in 1:30)(
  abline(a=fdakma_example$shift[i],b=fdakma_example$dilation[i], 
         col=fdakma_example_noalign_0der$labels[i])
)  





### COMPARE NUMBER OF CLUSTERS AND WARPING

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
