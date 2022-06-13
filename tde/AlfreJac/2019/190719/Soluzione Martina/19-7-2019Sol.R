# EXAM 19-7-2019


####################################
########### Exercise 1 #############
####################################

pines <- read.table('pines.txt', header=T)

# a) Perform a statistical test of level α = 0.01 to verify if the centre of the Pineta Dannunziana 
# can be assumed to be located in position Long = 14.2350, Lat = 42.4520. Report the p-value of 
# the test and verify its assumptions.

# H0: mu == c(14.2350,42.4520) vs H1:H0^C

n <- dim(pines)[1]
p <- dim(pines)[2]
alpha <- 0.01
mu0 <- c(14.2350,42.4520)

# verify the gaussianity assumption on my data:
mcshapiro.test(pines) # pvalue=0.097 -> ok

x.mean <- sapply(pines,mean)
x.cov <- cov(pines)
x.invcov <- solve(x.cov)

x.T2       <- n * (x.mean-mu0) %*% x.invcov %*% (x.mean-mu0) 
# Radius of the ellipsoid
cfr.fisher <- ((n-1)*p/(n-p))*qf(1-alpha,p,n-p)  # qf-> compute the quantile of F distribution 
# Test: 
x.T2 < cfr.fisher   # TRUE -> no statistical evidence to reject H0 at level alpha

# Compute the p-value 
P <- 1-pf(x.T2*(n-p)/((n-1)*p), p, n-p)
P # = 0.015

# So I can't reject H0.

# b) Estimate an elliptical region A that contains 99% of the pines. Report: the analytical 
# expression of the region, its centre, the direction and the length of the principal axes of the 
# ellipse. Report a qualitative plot of the region.

eigen(x.cov)

# Direction of the axes:
eigen(x.cov)$vectors

# Centre:
M <- mu0
M

# Radius of the ellipse:
r <- sqrt(qchisq(0.99,2))
r

# Length of the semi-axes:
r*sqrt(eigen(x.cov)$values)

# Plot
plot(pines, asp=1)
points(x.mean[1],x.mean[2],xlim=c(13,15), col='blue',pch=19,xlab='X.1',ylab='X.2',asp=1)
ellipse(M, x.cov, r, col = 'red', lty = 2, lwd=2)


####################################
########### Exercise 2 #############
####################################

buoys <- read.table('buoys.txt', header=T)

# a) Cluster the buoys based only on the GPS coordinates by using a hierarchical clustering method 
# (Euclidean distance and Ward linkage). Report a qualitative plot of the dendrogram and evaluate 
# the number of clusters you deem appropriate for the data. Report the GPS coordinates of the 
# centers of the clusters and their numerosity.

plot(buoys[,1:2]) # Note 3 different groups!

b.e <- dist(buoys[,1:2], method='euclidean')
b.ew <- hclust(b.e, method='ward.D')

# dendrogram:
plot(b.ew, main='euclidean-ward', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(b.ew, k=3)

cluster.ew <- cutree(b.ew, k=3) 
cluster.ew

# GPS coordinates of the centers of the clusters

m1 <- sapply(buoys[which(cluster.ew==1),1:2],mean)
#  Long      Lat 
# 14.68148 42.50352 
m2 <- sapply(buoys[which(cluster.ew==2),1:2],mean)
# Long      Lat 
# 14.40212 42.41839 
m3 <- sapply(buoys[which(cluster.ew==3),1:2],mean)
# Long      Lat 
# 14.40212 42.41839 

# numerosity of the clusters:
n1 <- length(which(cluster.ew==1)) # 348
n2 <- length(which(cluster.ew==2)) # 152
n3 <- length(which(cluster.ew==3)) # 252

# b) Assume the observations at different buoys to be independent. Perform a statistical test to 
# verify if there is a significant difference between the mean dissolved oxygen in the groups 
# identified at point (a). Formulate an appropriate model, state and verify the corresponding 
# assumptions. Comment the results.

O1 <- buoys[which(cluster.ew==1),3]
O2 <- buoys[which(cluster.ew==2),3]
O3 <- buoys[which(cluster.ew==3),3]


# H0: mu.O1 == mu.O2 == mu.O3 vs H1:H0^C
attach(buoys)
model <- aov(DO ~ cluster.ew)
summary(model)

# Assumption of the model
# 1) normality (univariate) in each group 
Ps <- c(shapiro.test(O1)$p,
        shapiro.test(O2)$p,
        shapiro.test(O3)$p)
Ps # ok

# 2) same covariance structure
Var <- c(var(O1),
         var(O2),
         var(O3)) 
Var

# There is no statistical evidence to say that there is difference between the mean dissolved oxygen 
# in the 3 groups.

# I can investigate which groups are more different than the other:
# Bonferroni confidence interval
g <- 3
n <- length(DO)
k <- g*(g-1)/2
alpha= 0.05

Mediag  <- tapply(DO, cluster.ew, mean)
SSres <- sum(residuals(model)^2)
S <- SSres/(n-g)
ng <- c(n1,n2,n3)

# CI for all the differences
ICrange=NULL
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    print(as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
    ICrange=rbind(ICrange,as.numeric(c(Mediag[i]-Mediag[j] - qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )),
                                       Mediag[i]-Mediag[j] + qt(1-alpha/(2*k), n-g) * sqrt( S * ( 1/ng[i] + 1/ng[j] )))))
  }}


# (1 - 2):  -0.02306509  0.02199431
# (1 - 3):  -0.06297574 -0.02464057
# (2 - 3):  -0.06707131 -0.01947422

# The first difference contains the value zero, so I can assess that there is more difference 
# between the first and the third groups and the second with the third one.


####################################
########### Exercise 3 #############
####################################

piadeina <- read.table('piadeina.txt', header=T)
attach(piadeina)

# a) Formulate a linear regression model for the total sales, as a function of all the other 
# variables. Include in the model a possible dependence of the total sales on the categorical 
# variable ‘day of the week’, but only in the intercept. Report the estimates of the 15 parameters
# of the model (the coefficients β0,...,β13 and the errors’ standard deviation σ). Analyse the 
# model residuals, and verify the assumptions of the model.

model <- lm(Sales ~ Bread.Sold + Wraps.Sold + Sandwich.Sold + Focaccia.Sold + Piadina.Sold 
            + Chips.Sold + Juices.Sold + Total.Soda.and.Coffee.Sold + Max.Daily.Temperature
            + Day.of.Week)
summary(model)

# Coefficients:
coeffs <- data.frame( beta0.0 = model$coefficients[1],
                      beta0.1 = model$coefficients[1] + model$coefficients[11],
                      beta0.2 = model$coefficients[1] + model$coefficients[12],
                      beta0.3 = model$coefficients[1] + model$coefficients[13],
                      beta0.4 = model$coefficients[1] + model$coefficients[14],
                      beta1 = model$coefficients[2],
                      beta2 = model$coefficients[3],
                      beta3 = model$coefficients[4],
                      beta4 = model$coefficients[5],
                      beta5 = model$coefficients[6],
                      beta6 = model$coefficients[7],
                      beta7 = model$coefficients[8],
                      beta8 = model$coefficients[9],
                      beta9 = model$coefficients[10],
                      sigma = sqrt(sum(residuals(model)^2)/model$df))
coeffs

# Test the gaussian hypothesis
shapiro.test(model$residuals) # pvalue= 0.12 -> ok
plot(model$residuals) # homoschedasticity

# b) Perform a variable selection through a Lasso method, by setting the parameter controlling 
# the penalization to λ = 5. Report the significant coefficients.

x <- model.matrix(Sales ~ Bread.Sold + Wraps.Sold + Sandwich.Sold + Focaccia.Sold + Piadina.Sold 
                  + Chips.Sold + Juices.Sold + Total.Soda.and.Coffee.Sold + Max.Daily.Temperature
                  + Day.of.Week)[,-1]
# Build the vector of response
y <- Sales

# Let's set a grid of candidate lambda's for the estimate
fit.lasso <- glmnet(x,y, lambda = 5)
summary(fit.lasso)
coef.lasso <- predict(fit.lasso, s=5, type = 'coefficients')[1:14,]
coef.lasso
# I can remove: bread, focaccia, juices, temperature, chips

#  remove these
x <- model.matrix(Sales ~ Wraps.Sold + Sandwich.Sold + Piadina.Sold 
                  + Total.Soda.and.Coffee.Sold
                  + Day.of.Week)[,-1]
# Build the vector of response
y <- Sales

# Let's set a grid of candidate lambda's for the estimate
fit.lasso <- glmnet(x,y, lambda = 5)
coef.lasso <- predict(fit.lasso, s=5, type = 'coefficients')[1:8,]
coef.lasso

# c) Optimize the parameter λ within the range [0,100] via cross-validation. Report the optimal
# λ and the corresponding estimated coefficients.

lambda.grid <- seq(0,100,length=100)
cv.lasso <- cv.glmnet(x,y,lambda=lambda.grid) # default: 10-fold CV

bestlam.lasso <- cv.lasso$lambda.min
bestlam.lasso # 2.020202

# Coefficients
coef.lasso.cv <- predict(cv.lasso, s=5, type = 'coefficients')[1:8,]
coef.lasso.cv


####################################
########### Exercise 4 #############
####################################

temp <- read.table('watertemp.txt', header=T)

# a) Perform a smoothing of the data through a projection over a Fourier basis with 45 basis 
# elements. Report the first 3 Fourier coefficients obtained at the Stations 1 and 2.
time <- 1:365

basis <- create.fourier.basis(rangeval=c(0,365),nbasis=45)
temp.fd <- Data2fd(y = t(as.matrix(temp[,1:365])),argvals = time,basisobj = basis)
plot.fd(temp.fd)

# Coefficients:
temp.fd$coefs[1:3,1:2]

# b) Perform a functional principal component analysis of the smoothed data obtained at point 
# (a). Report the variance explained along the first 5 functional principal components, a 
# qualitative plot of the first 3 eigenfunctions and the screeplot. Interpret the principal 
# components.

pca <- pca.fd(temp.fd,nharm=5,centerfns=TRUE)
names(pca)

# variance explained along the first 5 pc:
cumsum(pca$values)[1:5]/sum(pca$values)

# first two FPCs
layout(cbind(1,3))
plot(pca$harmonics[1,],col=1,ylab='FPC1',ylim=c(0,0.1))
plot(pca$harmonics[2,],col=2,ylab='FPC2',ylim=c(-0.1,0.1))
plot(pca$harmonics[3,],col=1,ylab='FPC3',ylim=c(-0.1,0.1))

# screeplot
plot(pca$values,xlab='j',ylab='Eigenvalues')
plot(cumsum(pca$values)/sum(pca$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

# Try to interprete the PCs:
# plot of the FPCs as perturbation of the mean
media <- mean.fd(temp.fd)

# first PC
plot(media,lwd=2,ylim=c(10,30),ylab='temperature',main='FPC1')
lines(media+pca$harmonics[1,]*sqrt(pca$values[1]), col=2)
lines(media-pca$harmonics[1,]*sqrt(pca$values[1]), col=3)
# variation in amplitude (more in the second period of the year than in the first)

# second PC -> contrast during the first and the second period of the year
plot(media,lwd=2,ylim=c(10,35),ylab='temperature',main='FPC2')
lines(media+pca$harmonics[2,]*sqrt(pca$values[2]), col=2)
lines(media-pca$harmonics[2,]*sqrt(pca$values[2]), col=3)

# third PC
plot(media,lwd=2,ylim=c(10,35),ylab='temperature',main='FPC3')
lines(media+pca$harmonics[3,]*sqrt(pca$values[3]), col=2)
lines(media-pca$harmonics[3,]*sqrt(pca$values[3]), col=3)
# So similar -> nothing important to explain

# c) Having reported a qualitative plot of the scores along the first 2 functional principal 
# components, use the categorical variable zone to further enhance the interpretations.

plot(pca$scores[,1],pca$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2) # 3 clusters!

zone <- as.factor(temp$Zone)
plot(pca$scores[,1],pca$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2, col = c('black', 'red', 'green')[zone])

# I can see that the different type of zone impact most on the first principal component, in fact 
# in the black group (deep zone) the scores are negative (decrease of temperature) and in the
# green one (surface zone) the scores are positive (increase of temperature). In the red one 
# (medium zone) we both have the effects. 
# Instead, there is no important effect of the variable zone in the second principal component.

# d) Propose a possible dimensionality reduction for the data and discuss the results.

# I can use the first principal component as dimensionality reduction, because it explained 
# the most percentage of variability.






