##### UTILITIES ESAME

rm(list=ls())
graphics.off()
load('mcshapiro.test.RData')
tilda = '~'

library(MASS) #qda
library(car) #vif
library(glmnet) #lasso
library(fda)
library(KernSmooth) #fda - smoothing
library(fdakma) #kmean alignment
library(heplots) #boxM
library(class) #knn

library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)  

#order -> black, red, green