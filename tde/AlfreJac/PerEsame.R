### COSE UTILI PER ESAME ###
## APPLIED STATISTICS ##

## Argomenti Lab ####

#### LAB 1
### Basic commands (scalars, vectors, matrices, and operations)
### Import/Export of dataframes
### Examples of univariate statistical analyses with plots
### Visualization of Multivariate Data
### Visualization of Categorical Data
### 3d plots, functions, "for" cycles
### save plots

#### LAB 2
### Probability density functions, Cumulative distribution functions, Quantiles
### Random number generation
### QQplots

#### LAB 3
### Principal Component Analysis

#### LAB 4
### Testing for multivariate normality

#### LAB 5
### Box-Cox transformations
### Tests and confidence regions for the mean of a multivariate Gaussian

#### LAB 6
### Test for the mean of paired multivariate Gaussian observations
### Test for repeated measures
### Test for two independent Gaussian populations

#### LAB 7
### One-way ANOVA and MANOVA
### Two-ways ANOVA and MANOVA

#### LAB 8
### Linear and Quadratic Discriminant Analysis
### Support Vector Machines

#### LAB 9
### Hierarchical clustering
### K-means clustering
### Exercises
### Multidimensional Scaling

#### LAB 10-11
### Linear models

## Spatial statistics

## Functional data analysis


## Libraries ####
library(mvtnorm)
library(mvnormtest)
library(rgl)
library(car)
library(MASS)
library(class) #clustering
library(e1071) #support vector machine
# Spatial statistics
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat) 
# Functional data analysis
library(fda)
library(KernSmooth)
library(fields)
library(fdakma)
load("/Users/martina/Downloads/Applied statistics/Labs/LAB_5/mcshapiro.test.RData")


## Aprire file .txt ####
nome <- read.table('nome.txt', header=T)

## 


