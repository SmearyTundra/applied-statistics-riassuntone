# EXAM 13/09/2018 ex 1

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

# The file IAMG.txt collects the data regarding the participation to the annual
# meetings of the International Association for Mathematical Geosciences (IAMG),
# in the last 50 years. For each meeting, it reports the number of registered
# participants, the number of oral presentations and the number of no-show
# (i.e., the number of registered participants that did not show up).
# Call X the vector whose components are the number of registered participants
# (X1), of oral presentations (X2) and of no-show (X3), and assume each meeting
# to be independent of the others.

IAMG <- read.table('IAMG.txt', header=T)
head(IAMG)

X <- IAMG
X1 <- IAMG[,1] # registered participants
X2 <- IAMG[,2] # oral presentations
X3 <- IAMG[,3] # no-show

# a) Build a confidence region (level 95%) for the mean of X. Characterize the
#    region by reporting its expression, its center, the direction of the axes
#    and the length of the semi-axes.

m <- sapply(X, mean)
S <- cov(X)
S.inv <- solve(S)

n <- dim(X)[1]
p <- dim(X)[2]

alpha <- 0.05
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha, p, n-p)

# Analytical expression:
# El(x) = {x in R^p s.t. (x-M)'*S.inv*(x-M)<=r^2}

# Center:
m
# Registered       Talk    No.show 
#     249.44     119.04      23.12 

# Direction of the principal axes:
eigen(S/n)$vectors
# 0.98711610  0.1554660 -0.03784352
# 0.09129046 -0.7414562 -0.66476221
# 0.13140723 -0.6527427  0.74609588

# Length of the semi-axes:
r <- sqrt(cfr.fisher)
r*sqrt(eigen(S/n)$values)
# 9.914474 4.218786 2.708017

# Gaussianity assumption
mcshapiro.test(X) # pvalue=0.228 -> OK

# ____________________________________________________________________________

# b) Build three T2-simultaneous confidence intervals (level 95%) for: the mean
#    number of registered participants, the mean number of oral presentations
#    and the mean number of no-show.

IC.T2 <- cbind(m - sqrt(diag(S)/n*cfr.fisher),
               m,
               m + sqrt(diag(S)/n*cfr.fisher))
IC.T2
# Registered 239.63077 249.44 259.24923
# Talk       115.31917 119.04 122.76083
# No.show     19.46447  23.12  26.77553

# ___________________________________________________________________________

# c) Perform a test of level 95% to verify the hypothesis according to which,
#    in mean, only 90% of the registered participants actually show up at IAMG
#    meetings.

# Test:
# H0: 0.10*mu.Registered==mu.No.show   vs   H1: 0.10*mu.Registered!=mu.No.show

# Test on a linear combination of the means

a <- c(0.10, 0, -1)
delta.0 <- 0

Ma <- t(a)%*%m
Sa <- t(a)%*%S%*%a
Sa.inv <- solve(Sa)

T2 <- n*(Ma-delta.0)%*%Sa.inv%*%(Ma-delta.0)
T2 # 2.394583

T2 > cfr.fisher # FALSE -> Do not reject H0 -> We can confirm that only 90% of
# participants actually show up at IAMG

