# EXAM 15/06/2020 ex 1

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

library(car)

# The file pollution.txt reports the measurements of two air pollutants (PM2.5
# and PM10) collected in 100 days by an air-quality monitoring site of a Chinese
# city.

pollution <- read.table('pollution.txt', header=T)
head(pollution)

# a) Perform a statistical test at level 5% to verify if the mean of the two air
#    pollutants is significantly different from (50,50). Verify the assumptions
#    required to perform the test and provide the p-value.

# TEST FOR THE MEAN OF PAIRED MULTIVARIATE GAUSSIAN OBSERVATIONS

# Test: H0: mu==delta.0   vs   H1: mu!=delta.0

# Check Gaussianity
mcshapiro.test(pollution) # 0.3008

alpha <- 0.05
delta.0 <- c(50,50)

n <- dim(pollution)[1]
p <- dim(pollution)[2]

M <- sapply(pollution, mean)
S <- cov(pollution)
S.inv <- solve(S)

T2 <- n*(M-delta.0)%*%S.inv%*%(M-delta.0)
T2 # 262.631
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha, p, n-p)
cfr.fisher # 6.241451
T2 > cfr.fisher # TRUE -> Reject H0 -> There is statistical evidence that the
# mean of the two air pollutants is significantly different from (50,50)

# _________________________________________________________________________

# b) Find an elliptical confidence region at level 95% for the mean value of the
#    two air pollutants. Provide a plot and report: the analytical expression of
#    the region, its centre, the direction and the length of the principal axes
#    of the ellipse.

# Analytical expression:
# El(x) = {x in R2 s.t. (x-M)'S.inv(x-M)<=r^2}

# Centre:
M
# PM2.5     PM10 
# 164.9005  150.4848 

# Directions on the principal axes:
eigen(S)$vectors
# 0.6291343 -0.7772966
# 0.7772966  0.6291343

# Length of the principal axes:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(S)$values)
# 262.95269  85.25918

# Plot:
plot(pollution, pch=19, xlim=c(-20, 350), ylim=c(-50, 400))
ellipse(center = M, shape = S/n, radius = r) # ???

# ________________________________________________________________________

# c) Does the vector (50,50) lies within or outside the region identified at
#    point(b)? Comment the answer highlighting the connection to the conclusion
#    of point (a).

points(50,50) 
# It lies outside of the elliptical confidence region.
# It makes sense, since in point (a) we rejected H0 at 95%, i.e. we rejected the
# possibility to state that point (50,50) could be equal to the sample mean,
# which is the center of our ellipse.

# __________________________________________________________________________

# d) Provide two T2 simultaneous confidence intervals (global confidence 95%)
#    for the mean of PM2.5 and the mean of PM10.

IC.T2 <- cbind(M - sqrt(diag(S)/n*cfr.fisher),
            M,
            M + sqrt(diag(S)/n*cfr.fisher))
IC.T2
# PM2.5 147.0792 164.9005 182.7218
# PM10  129.3534 150.4848 171.6161