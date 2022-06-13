
# EXAM 14/9/2017

####################################
########## EXERCISE 1 ##############
####################################

sunchair <- read.table('sunchair.txt', header=T)

# a) Through an appropriate statistical test, verify if there exist a significant variation 
# in the mean prices of a sun chair during the year. Introduce and verify the needed assumptions.

# Use repeated measures framework

# gaussian hypothesis:
mcshapiro.test(sunchair) # 0.1916 -> ok

attach(sunchair)

n <- dim(sunchair)[1]
g <- dim(sunchair)[2]
n_tot <- g*n

labels <- as.factor(c(rep("Mar.May", n), rep("Jun.July", n), rep("Aug.Oct", n), rep("Nov.Feb", n)))
new.data <- data.frame(prices = c(Mar.May, Jun.July,Aug.Oct,Nov.Feb),
                       period = labels)
attach(new.data)
fit.aov <- aov(prices ~ period)
summary(fitaov)
# there is statistical evidence to assess that there is difference!

# b) Use four Bonferroni intervals (global level 95%) to describe the dynamic of the mean price 
# of a sun chair. Based on your analyses, suggest the best period to buy a sun chair and its 
# expected price.

# Construct Bonferroni interval for every period, and not for the increments -> I think is better
# to understand which period is the best to buy a sunchair
k <- 4  
cfr.t <- qt(1-alpha/(2*k),n-1)
M <- sapply(sunchair,mean)
S <- cov(sunchair)

IC.BF <- cbind( M - cfr.t*sqrt(diag(S)/n) , M, M + cfr.t*sqrt(diag(S)/n) )
IC.BF

#             inf     center      sup
# Mar.May  41.54536 42.60177 43.65819
# Jun.July 49.63851 50.44403 51.24956
# Aug.Oct  44.11548 44.69871 45.28194
# Nov.Feb  37.69024 38.66500 39.63976

# The best period to buy a sunchair is November-February, and the expected prices is 39 euros


####################################
########## EXERCISE 2 ##############
####################################

olives <- read.table('olives.txt', header=T)
attach(olives)

# Divide the two group
olives1 <- olives[which(Restaurant=='Dalla Luigina'),]
olives2 <- olives[-which(Restaurant=='Dalla Luigina'),]

# a) Is there a significant difference (level 95%) in the mean of the total weight, 
# the filling weight and the breading weight of the olives served in the two restaurants? 
# Introduce and verify the needed assumptions.

# we build the data
t1 <- olives1[,-4]
t2 <- olives2[,-4]

# Assumption: two independent gaussian population -> OK
mcshapiro.test(t1) # pvalue = 0.9724
mcshapiro.test(t2) # pvalue = 0.1964

n1 <- dim(t1)[1] # n1=40
n2 <- dim(t2)[1] # n2=36
p  <- dim(t1)[2] # p=3

# compute the sample mean, covariance matrices and the matrix Spooled

t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)

# Test H0: mu1 == mu2  vs  H1: mu1 != mu2
# i.e.,
# Test H0: mu1-mu2 == c(0,0)  vs  H1: mu1-mu2 != c(0,0)

alpha   <- .05
delta.0 <- c(0,0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (t1.mean-t2.mean-delta.0) %*% Spinv %*% (t1.mean-t2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # FALSE: reject H0

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  # = 0 -> reject H0: there is difference in mean

# b) Provide T2 intervals for the mean difference between the total weight, the filling weight 
# and the breading weight of the olives served Dalla Luigina and at Caffe Muletti. Comment the results.

IC.T2 <- data.frame(inf = t1.mean-t2.mean-sqrt(cfr.fisher*diag(Sp)*(1/n1+1/n2)),center = t1.mean - t2.mean, sup =t1.mean-t2.mean+sqrt(cfr.fisher*diag(Sp)*(1/n1+1/n2))) 
IC.T2

# 1:dalla luigina, 2:muletti -> All negative intervals: Muletti makes bigger olives! In particular,
# the filling weigth is bigger


####################################
########## EXERCISE 3 ##############
####################################

knossos <- read.table('knossos.txt', header=T)

# a) Identify two clusters of locations through a hierarchical clustering method (Euclidean 
# distance and complete linkage). Report the estimate of the mean within the groups, their size, 
# and compute the cophenetic coefficient

knossos.e <- dist(knossos, method='euclidean')
knossos.ec <- hclust(knossos.e, method='complete')
cluster.ec <- cutree(knossos.ec, k=2) 
cluster.ec

knossos1 <- knossos[which(cluster.ec==1),]
knossos2 <- knossos[which(cluster.ec==2),]

# mean within the two groups
mean1 <- sapply(knossos1, mean) 
#   X        Y 
# 2.962927 1.503984 

mean2 <- sapply(knossos2, mean)
#   X           Y 
# 0.03269231 -0.02320513 

# size of the groups
n1 <- length(knossos1[,1])
n2 <- length(knossos1[,2])

# cophenetic coefficient
coph.ec <- cophenetic(knossos.ec)
ec <- cor(knossos.e, coph.ec) # 0.8709168

# b) Assume the identified groups to be independent. Having introduced and verified the 
# needed assumptions, test the hypothesis according to which only one archeological site exists 
# in the Knossos area. Write a report of max 3 lines to the archeologists summarising the results 
# of the analysis.

plot(knossos.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(knossos.ec, k=2)

# H0: exists only 1 cluster vs H1: H0^C
# i.e., the two groups have no difference in mean!
# H0: mean.knossos1 - mean.knossos2 == 0 vs H1:H0^C

# gaussian assumptions
mcshapiro.test(knossos1) # pvalue=0.62 -> ok
mcshapiro.test(knossos2) # pvalue=0.58 -> ok

p <- dim(knossos1)[2]
k1.mean <- sapply(knossos1,mean)
k2.mean <- sapply(knossos2,mean)
k1.cov  <-  cov(knossos1)
k2.cov  <-  cov(knossos2)
Sp      <- ((n1-1)*k1.cov + (n2-1)*k2.cov)/(n1+n2-2)

# we compare the matrices
list(S1=k1.cov, S2=k2.cov, Spooled=Sp) # same covariance structures

alpha   <- .05
delta.0 <- c(0,0)
Spinv   <- solve(Sp)

T2 <- n1*n2/(n1+n2) * (k1.mean-k2.mean-delta.0) %*% Spinv %*% (k1.mean-k2.mean-delta.0)

cfr.fisher <- (p*(n1+n2-2)/(n1+n2-1-p))*qf(1-alpha,p,n1+n2-1-p)
T2 < cfr.fisher # FALSE: reject H0 at level 5% -> there is difference!

P <- 1 - pf(T2/(p*(n1+n2-2)/(n1+n2-1-p)), p, n1+n2-1-p)
P  # reject H0

# So, there exist two different groups of archeological sites according to my test.


####################################
########## EXERCISE 4 ##############
####################################

tide <- read.table('tide.txt', header=T)
attach(tide)

# a) Estimate the five parameters of the model. Report the estimates of βi, i ∈ {1,2,3}, and σ.

# Make the model
x1 <- I(sin(2*pi/28*t))
x2 <- I(sin(pi/365*(t-82)))
x3 <- t

model <- lm(h ~ x1 +x2 + x3)
summary(model)
shapiro.test(model$residuals) # pvalue=0.24 -> ok

coefficients(model)
# (Intercept)          x1          x2          x3 
# 66.49527511 19.21718913  1.69733063  0.02371382 

sqrt(sum(residuals(model)^2)/model$df) # sigma
# 8.904655

# b) Having introduced and verified the appropriate assumptions, perform two statistical tests to verify if
# - the mean sea level is influenced by the periodic components;
# - the mean sea level depends on the global increase of the sea level.

linearHypothesis(model, c(0,0,1,0), 0)
# not influenced by periodic components

new.model <-  lm(h ~ x1 + x3)
summary(new.model)

linearHypothesis(new.model, c(0,0,1), 0)
# influenced 

# c) Based on point (b), propose a reduced model and estimate its parameters.

# Reduce model
model_reduce <- lm(h ~ x1 + x3)
summary(model_reduce)

# test gaussian assumption on the residuals
shapiro.test(model_reduce$residuals) # pvalue= 0.18 -> ok

coefficients(model_reduce)
# (Intercept)          x1 
# 71.80895    18.96390 

sqrt(sum(residuals(model_reduce)^2)/model_reduce$df)
# 9.45829

# d) Based on model (c), provide two prediction intervals (global level 90%) for the sea level at 
# 17:00 of 20th September 2017 (day 263 of 2017) and of 1st December 2017 (day 335 of 2017). 
# Comment the results knowing that in Venice high-water is expected whenever the sea level is higher 
# than 90 cm

Z01.new <- data.frame(x1 = I(sin(2*pi/28*263)), x3 = 263) #new dataframe

Pred1 <- predict(model_reduce, Z01.new, interval='prediction', level=1-0.1)  
Pred1

# fit      lwr      upr
# 86.08586 71.30669 100.865

Z02.new <- data.frame(x1= I(sin(2*pi/28*335)), x3 = 335)

Pred2 <- predict(model_reduce, Z02.new, interval='prediction', level=1-0.1)  
Pred2

# fit      lwr      upr
# 72.0717 57.26924 86.87416

# the first prediction interval (at day 263) containes 90 cm level (Venice high-water!), 
# the second not
