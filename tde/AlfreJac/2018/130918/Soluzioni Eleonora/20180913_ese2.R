# EXAM 13/09/2018 ex 2

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

# The file Waiting.txt reports the waiting times for food - i.e., the time
# between the order of a course and the service - in 180 restaurants in Romania.
# The dataset also reports the type of course (starter, main course or dessert)
# and the location of the restaurant (Iasi or Bucarest).

data <- read.table('Waiting.txt', header=T)
head(data)
attach(data)

# a) Propose a complete ANOVA model for the waiting time as a function of the
#    factors course (starters, main course or dessert) and city (Iasi or
#    Bucarest). Report and verify the assumptions of the model.

# Two-ways ANOVA

# Variable: waiting time
# Factor 1: course (g=3)
# Factor 2: city (b=2)

# Model:
# X.ijk = mu + tau.i + beta.j + gamma.ij + eps.ijk,   eps ~ N(0, sigma^2)
# i=1,2,3 -> Factor course
# j=1,2 -> Factor city

g <- 3
b <- 2
n <- dim(data)[1]/(g*b)
N <- n*g*b

# Assumptions: 
# 1) Gaussianity (univariate) in each group
shapiro.test(waiting[course=='Starter' & city=='Iasi']) # pvalue = 0.7913 -> OK
shapiro.test(waiting[course=='Main' & city=='Iasi']) # pvalue = 0.4336 -> OK 
shapiro.test(waiting[course=='Dessert' & city=='Iasi']) # pvalue = 0.4099 - > OK
shapiro.test(waiting[course=='Starter' & city=='Bucarest']) # pvalue = 0.152 -> OK
shapiro.test(waiting[course=='Main' & city=='Bucarest']) # pvalue = 0.2473 -> OK
shapiro.test(waiting[course=='Dessert' & city=='Bucarest']) # pvalue = 0.8573 -> OK
# 2) Same covariance structure (homogeneity of variances)
bartlett.test(waiting, course) # pvalue = 0.1509 -> OK
bartlett.test(waiting, city) # pvalue = 0.647 -> OK

fit <- aov(waiting ~ course + city + course:city)

# ___________________________________________________________________________

# b) Comment on the significance of the factors and of their interaction.
#    If needed, propose a reduced model.

summary.aov(fit)
# In the one-at-a-time tests reported in the summary, we have:
# Factor course -> Significant (very small pvalue)
# Factor city -> Not significant (pvalue=0.907 -> very big)
# Interaction term -> Not significant (pvalue=0.573 -> very big)
# -> So we can try to reduce the model

# First we remove interaction:

# Model: Two-ways ANOVA (additive)
# X.ijk = mu + tau.i + beta.j + eps.ijk,   eps ~ N(0, sigma^2)
# i=1,2,3 -> Factor course
# j=1,2 -> Factor city
fit2 <- aov(waiting ~ course + city)
summary.aov(fit2)
# Factor city still has a very big pvalue

# We remove factor city:

# Model: One-way ANOVA
# X.ijk = mu + tau.i + eps.ik,   eps ~ N(0, sigma^2)
# i=1,2,3 -> Factor course
fit3 <- aov(waiting ~ course)
summary.aov(fit3)

# __________________________________________________________________________

# c) Build Bonferroni confidence intervals (global level 95%) for the mean
#    differences between the waiting times in the groups identified at point (b),
#    and for the variances of the waiting times within the groups. Comment the
#    results.

alpha <- 0.05
k <- 4 # 3 CI for the differences of the means and 1 CI for the variances
nn <- 60

i1 <- which(course=='Starter')
i2 <- which(course=='Main')
i3 <- which(course=='Dessert')

n1 <- length(i1)
n2 <- length(i2)
n3 <- length(i3)


alpha <- 0.05
k <- g*(g-1)/2 # 3 ci for the means


DF <- fit3$df.residual
qT <- qt(1-alpha/(2*k), DF)

SSres <- sum(residuals(fit3)^2)
S <- SSres/DF
#m  <- sapply(data, mean)         # estimates mu
m1 <- mean(data[i1, 1])    # estimates mu.1=mu+tau.1
m2 <- mean(data[i2, 1])    # estimates mu.2=mu+tau.2
m3 <- mean(data[i3, 1])    # estimates mu.3=mu+tau.3

inf12 <- m1-m2 - qT * sqrt(S*(1/n1+1/n2) )
sup12 <- m1-m2 + qT * sqrt(S*(1/n1+1/n2) )
inf13 <- m1-m3 - qT * sqrt(S*(1/n1+1/n3) )
sup13 <- m1-m3 + qT * sqrt(S*(1/n1+1/n3) )
inf23 <- m2-m3 - qT * sqrt(S*(1/n2+1/n3) )
sup23 <- m2-m3 + qT * sqrt(S*(1/n2+1/n3) )

IC.BF <- list(Starter.Main=cbind(inf12, sup12),
              Starter.Dessert=cbind(inf13, sup13),
              Main.Dessert=cbind(inf23, sup23))
IC.BF
# $Starter.Main        $Starter.Dessert    $Main.Dessert
# inf12     sup12      inf13    sup13      inf23    sup23
# -5.512386 1.479053   45.47095 52.46239   47.48761 54.47905

# The difference in the mean waiting time is significant in Starter vs Dessert
# and in Main vs Dessert (both positive intervals), while it's not in
# Starter vs Main (in this case the confidence interval contains 0). 
# This means that Dessert have a mean waiting time which is a lot lower than the
# other courses.
