# EXAM 14/9/2017 ex 2

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

# The province of Ascoli Piceno is well-known for the olive all'ascolana
# (i.e., breaded and fried stuffed olives), hereafter named just olives.
# The file olives.txt reports the total weight [g], the filling weight [g] and
# the weight of the fried breading of 40 olives served in the restaurant
# Dalla Luigina, and of 36 olives served at the Caffe Muletti, in Ascoli Piceno.

olives <- read.table('olives.txt', header=T)
head(olives)

# a) Is there a significant difference (level 95%) in the mean of the total
#    weight, the filling weight and the breading weight of the olives served
#    in the two restaurants? Introduce and verify the needed assumptions.

# Divide the 2 groups:
olives1 <- olives[which(olives$Restaurant=='Dalla Luigina'),]
olives2 <- olives[-which(olives$Restaurant=='Dalla Luigina'),]
olives1 <- olives1[,-4]
olives2 <- olives2[,-4]

# TEST FOR THE MEANS OF 2 INDEPENDENT GAUSSIAN POPULATIONS

# Check the Gaussianity assumption in each group:
mcshapiro.test(olives1) # pvalue=0.976 -> OK
mcshapiro.test(olives2) # pvalue=0.2136 -> OK

n1 <- dim(olives1)[1]
n2 <- dim(olives2)[1]
p <- dim(olives1)[2]

m1 <- sapply(olives1, mean)
m2 <- sapply(olives2, mean)
S1 <- cov(olives1)
S2 <- cov(olives2)
Sp <- ((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
Sp.inv <- solve(Sp)

# Test: H0: mu.1==mu.2   vs   H1: mu.1!=mu.2
alpha <- 0.05
delta.0 <- c(0,0,0)

T2 <- n1*n2/(n1+n2)*(m1-m2-delta.0)%*%Sp.inv%*%(m1-m2-delta.0)
T2 # 532.7752
cfr.fisher <- (n1+n2-2)*p/(n1+n2-1-p)*qf(1-alpha, p, n1+n2-1-p)
cfr.fisher # 8.423072

T2 > cfr.fisher # TRUE -> Reject H0 -> There is a difference in mean between
# the two restaurants

# _______________________________________________________________________

# b) Provide T2 intervals for the mean difference between the total weight,
#    the filling weight and the breading weight of the olives served
#    Dalla Luigina and at Caffe Muletti. Comment the results.

IC.T2 <- cbind(m1-m2-delta.0 - sqrt(cfr.fisher*diag(Sp)*(1/n1+1/n2)),
               m1-m2-delta.0,
               m1-m2-delta.0 + sqrt(cfr.fisher*diag(Sp)*(1/n1+1/n2)))
IC.T2
# Tot   -9.9636729 -8.7677778 -7.5718826
# Fill  -6.6147750 -5.6411389 -4.6675027
# Bread -0.9277802 -0.6386667 -0.3495531

# All negative intervals -> Caffe Muletti makes bigger olives: breading in Caffe
# Muletti weights a little more than in Dalla Luigina, but the most relevant
# difference is observable in the filling, which is much bigger.