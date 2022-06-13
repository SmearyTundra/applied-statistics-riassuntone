## EXAM 3-7-2017

####################################
########### EXERCISE 1 #############
####################################

kimono <- read.table('kimono.txt', header=T)
attach(kimono)

# a) Formulate an ANOVA model for the value of a kimono as a function of the factors city 
# (Kyoto or Tokyo) and type (tailor-made, ready-to-wear). Verify the assumptions of the model.

# two-ways ANOVA
fit.aov <- aov(value ~ city + type)
summary.aov(fit.aov)

# Verify the gaussianity assumptions -> OK
Ps <- c(shapiro.test(value[ city=='Kyoto' & type=='hand-made' ])$p,
        shapiro.test(value[ city=='Tokyo' & type=='hand-made' ])$p,
        shapiro.test(value[ city=='Kyoto' & type=='ready-to-use' ])$p,
        shapiro.test(value[ city=='Tokyo' & type=='ready-to-use' ])$p)
Ps
# homogeneity of variances -> OK
bartlett.test(list(value[ city=='Kyoto' & type=='hand-made' ],
                   value[ city=='Tokyo' & type=='hand-made' ],
                   value[ city=='Kyoto' & type=='ready-to-use' ],
                   value[ city=='Tokyo' & type=='ready-to-use' ]))

# b) Through appropriate statistical tests, propose a reduced model
summary(fit.aov)
# looking at the p-value of the one-at-the-time tests on the regressors, I can assess that the variable "city"
# is not significant for my model

fit.aov.reduce <- aov(value ~ type)
summary(fit.aov.reduce)

# c) Provide Bonferroni intervals (global level 95%) for the differences between the mean value of kimonos 
# belonging to the homogeneous groups identified by the model at point (b).

hand.made <- value[type=='hand-made']
ready.to.use <- value[type=='ready-to-use']
n1 <- length(hand.made)
n2 <- length(ready.to.use)
# n1 = n2 = n -> 264
alpha <- .05

n <- dim(kimono)[1]/2
g <- 2
k <- g*(g-1)/2 # 1

mean.val.type  <- tapply(value, type, mean) 

SSres <- sum(residuals(fit.aov.reduce)^2)

IC <- c(mean.val.type[1] - mean.val.type[2] - qt(1-alpha/(2*k), n*g-1) * sqrt(SSres/(n*g-1) * (1/(n) + 1/(n))), 
        mean.val.type[1] - mean.val.type[2] + qt(1-alpha/(2*k), n*g-1) * sqrt(SSres/(n*g-1) * (1/(n) + 1/(n))))
names(IC) <- c('Inf', 'Sup')
IC 
# Inf       Sup 
# -15.18700 -14.63194 

detach(kimono)


####################################
############ EXERCISE 2 ############
####################################

# a) Perform a statistical test to verify if there is evidence of an impact of Hanami on the
# mean amount of rice, sashimi, vegetables and okashi in families bent ̄o’s. Verify the needed assumptions.

bento <- read.table('bento.txt', header=T)
attach(bento)

D <- data.frame(D.rice = rice_hanami - rice_nohanami, 
                D.sashimi = sashimi_hanami- sashimi_nohanami, 
                D.veg= vegetables_hanami- vegetables_nohanami, 
                D.okashi= okashi_hanami - okashi_nohanami)

### T2 Hotelling Test 
# H0: D.mean == c(0,0,0,0) vs H1: D.mean != c(0,0,0,0)

# Test the Gaussian assumption -> p=0.9 -> OK
mcshapiro.test(D)

n <- dim(D)[1]  # 32
p <- dim(D)[2]  #  4

D.mean   <- sapply(D,mean)
D.cov    <- cov(D)
D.invcov <- solve(D.cov)

alpha   <- .05
delta.0 <- c(0,0,0,0)

#test statistic -> Mahalanobis distance between sample mean and delta0
D.T2 <- n * (D.mean-delta.0) %*% D.invcov %*% (D.mean-delta.0)
D.T2

cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)
cfr.fisher

D.T2 < cfr.fisher # FALSE: we reject H0 at level 5%

# we compute the p-value
P <- 1-pf(D.T2*(n-p)/(p*(n-1)), p, n-p)
P
# reject H0 at 5% 
# -> So there is no statistical evidence to say that the mean value of hanomi and nohanomi is equal

# b) Provide four T2 simultaneous confidence intervals (global confidence 95%) for the increase 
# in the mean consumption of rice, sashimi, vegetables and okashi in correspondence of the bloom 
# of cherry blossoms. Comment the results.

IC.T2. <- c( D.mean-sqrt(cfr.fisher*diag(D.cov)/n) , D.mean, D.mean+sqrt(cfr.fisher*diag(D.cov)/n) )
T2 <- rbind(IC.T2.rice, IC.T2.sashimi, IC.T2.veg, IC.T2.okashi)
dimnames(T2)[[2]] <- c('inf','center','sup')
T2

#                   inf     center       sup
# IC.T2.rice    -37.34390  -2.145937  33.05202
# IC.T2.sashimi 107.51324 123.577812 139.64238
# IC.T2.veg     -12.40493  11.463438  35.33181
# IC.T2.okashi  111.98999 118.150000 124.31001
#

## Comment: for rice and vegetables, the confidence intervals contain the zero value, so we can assess
# that in mean the difference can be zero. For the other not.


####################################
########## EXERCISE 3 ##############
####################################

geisha <- read.table('geisha.txt', header=T)

# a) Use a hierarchical clustering method based on Euclidean distance and single linkage 
# to identify two groups of data (i.e., successful and unsuccessful tours). Report the centers of the 
# clusters, the size of the clusters, the cophenetic coefficient and a qualitative plot of the results.

# euclidean distance
geisha.e <- dist(geisha, method='euclidean')

# clustering
geisha.es <- hclust(geisha.e, method='single')

cluster.es <- cutree(geisha.es, k=2)

# centers of the clusters
mean(which(cluster.es==1)) # 0.67296
mean(which(cluster.es==2)) # 53

# size of the clusters
length(which(cluster.es==1)) # 159
length(which(cluster.es==2)) # 1 -> 53

# cophenetic coefficient
coph.es <- cophenetic(geisha.es)
es <- cor(geisha.e, coph.es)  # 0.8781562

# qualitative plot of the results
plot(geisha, col=ifelse(cluster.es==1,'red','blue'), pch=19)
# I can see that there are another two possible clusters to be indentified!
# I'm not satisfied, I choose another procedure of clustering

# b) Evaluate the quality of the clustering at point (a) and, in case you deem it unsatisfactory, 
# repeat the procedure with another linkage at your choice.

# I use complete linkage:
geisha.ec <- hclust(geisha.e, method='complete')
cluster.ec <- cutree(geisha.ec, k=2)

plot(geisha, col=ifelse(cluster.ec==1,'red','blue'), pch=19)

# now I can see two different clusters! Satisfied

# c) Identify the successful tours with the smaller group found with the clustering method at 
# point (b). Having introduced and verified the needed assumptions, provide 4 Bonferroni 
# intervals (global level 90%) for the difference in the mean characteristics of successful 
# and unsuccessful tours, and for the mean characteristics of a successful tour.

length(which(cluster.ec==1)) # 75 -> successful cluster!
length(which(cluster.ec==2)) # 85

succ <- which(cluster.ec==1)

t1 <- geisha[succ,]
t2<- geisha[-succ,]

n1 <- dim(t1)[1] # 75
n2 <- dim(t2)[1] # 85
p  <- dim(t1)[2] # 2

# compute the sample mean, covariance matrices and the matrix Spooled
t1.mean <- sapply(t1,mean)
t2.mean <- sapply(t2,mean)
t1.cov  <-  cov(t1)
t2.cov  <-  cov(t2)
Sp      <- ((n1-1)*t1.cov + (n2-1)*t2.cov)/(n1+n2-2)

# assumption -> gaussianity
mcshapiro.test(t1) #ok
mcshapiro.test(t2) #ok

alpha <- 0.1
IC <- cbind(t2.mean-t1.mean - sqrt(diag(Sp)*(1/n1+1/n2)) * (p*(n1+n2-2)/(n1+n2-1-p))*qt(1 - alpha/(p*2), n1+n2-2),
            t2.mean-t1.mean,
            t2.mean-t1.mean + sqrt(diag(Sp)*(1/n1+1/n2)) * (p*(n1+n2-2)/(n1+n2-1-p))*qt(1 - alpha/(p*2), n1+n2-2))

#           inf       center    sup
# duration -47.90682 -45.20690 -42.50698
# time      13.64587  16.69882  19.75177

# intervals for the mean characteristics of a successful tour
IC <- cbind(t1.mean - sqrt(diag(t1.cov)*(1/n1)) * qt(1 - alpha/(p*2), n1-1),
            t1.mean,
            t1.mean + sqrt(diag(t1.cov)*(1/n1)) * qt(1 - alpha/(p*2), n1-1))

#           inf     center    sup
# duration 88.14106 90.23867 92.33627
# time     42.68666 45.20000 47.71334

# d) Comment the results at point (c) and suggest a successful strategy for Geisha hunting.

####################################
########## EXERCISE 4 ##############
####################################

garden <- read.table('garden.txt', header=T)
attach(garden)

# a) Estimate the 6 parameters of the model and verify the model assumptions. 
# Evaluate the residuals of the model.

fit <- lm(extension ~ carps + maple + cherry + stones)
summary(fit)

# Estimate the parameters
coefficients(fit)  # beta_i
residuals(fit)

# Assumptions on the residuals
shapiro.test(fit$residuals) # p-value = 0.8567 -> OK

# b) Perform two statistical tests to verify if
# - there is statistical evidence of a dependence of the mean garden extension on the number 
# of maple or cherry trees;
# - there is statistical evidence of a dependence of the mean garden extension on lake elements 
# (stones, carps).

linearHypothesis(fit, c(0,0,0,1,0), 0)  #test hypothesis on linear combination of vector beta
# Accept that cherry is not significant
new.model <- lm(extension ~ carps + maple + stones)
summary(new.model)

linearHypothesis(new.model, c(0,0,1,0), 0)
# Accept that maple is significant

linearHypothesis(new.model, rbind(c(0,1,0,0), c(0,0,0,1)), c(0,0))
# Reject null hypothesis -> the two cannot be take out

# c) Based on the results at point (b), comment on possible model weaknesses and, if needed, 
#reduce the dimensionality of the regressors. Comments the analysis and interpret the results.

# Possible weaknesses -> I reject the hypothesis that lake elements (together) are significant.
# But if I prove to remove one at the time:
linearHypothesis(new.model, c(0,1,0,0), 0)
# carps is not significant!
new.model <- lm(extension ~ maple + stones)

linearHypothesis(new.model, c(0,0,1), 0)
# stones is significant

# REDUCED MODEL
fit_reduce <- lm(extension ~ maple + stones)
summary(fit_reduce)
# now all the regressors are signfiicant

# d) Update the estimates of the parameters using the analysis at point (c).
coefficients(fit_reduce)
residuals(fit_reduce)

