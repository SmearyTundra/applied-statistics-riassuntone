# EXAM 14/9/2017 ex 1

setwd("D:/Dati/Ele/poli/anno4/applied statistics/OldExams/tde")
load("D:/Dati/Ele/poli/anno4/applied statistics/R lab/LAB_5/mcshapiro.test.RData")

# The file sunchair.txt collects the prices [e] of 62 sun chairs of comparable
# characteristics sold by a famous e-commerce company, in 4 different periods of
# the year: Mar-May, Jun-July, Aug-Oct, Nov-Feb.
# Assume the prices of different chairs to be independent.

sunchair <- read.table('sunchair.txt', header=T)
head(sunchair)

# a) Through an appropriate statistical test, verify if there exist a significant
#    variation in the mean prices of a sunchair during the year. Introduce and
#    verify the needed assumptions.

# REPEATED MEASURES FRAMEWORK

# Check Gaussian assumption:
mcshapiro.test(sunchair) # pvalue=0.1956 -> OK

# Test: H0: mu1=m12=mu3=mu4   vs   H1: Not constant mean price during the year
# i.e.  H0: C*mu==0   vs   H1: C*mu!=0

n <- dim(sunchair)[1]
q <- dim(sunchair)[2]

M <- sapply(sunchair, mean)
S <- cov(sunchair)

# Contrast matrix (successive increments)
C <- matrix(c(-1,1,0,0,
              0,-1,1,0,
              0,0,-1,1), 3, 4, byrow=T)

delta.0 <- c(0,0,0)
Md <- C%*%M
Sd <- C%*%S%*%t(C)
Sd.inv <- solve(Sd)

T2 <- n*t(Md-delta.0)%*%Sd.inv%*%(Md-delta.0)
P <- 1-pf(T2*(n-(q-1))/((q-1)*(n-1)), q-1, n-(q-1))
P
# Pvalue=0 -> Reject H0 for any alpha -> There is a significant variation of 
# prices during the year

# ___________________________________________________________________________

# b) Use four Bonferroni intervals (global level 95%) to describe the dynamic of
#    the mean price of a sun chair. Based on your analyses, suggest the best
#    period to buy a sun chair and its expected price.

alpha <- 0.05
k <-4
cfr.t <- qt(1-alpha/(2*k), n-1)

IC.BF <- cbind(M - cfr.t*sqrt(diag(S)/n),
               M,
               M + cfr.t*sqrt(diag(S)/n))
IC.BF
# Mar.May  41.49724 42.60177 43.70631
# Jun.July 49.60182 50.44403 51.28625
# Aug.Oct  44.08892 44.69871 45.30850
# Nov.Feb  37.64584 38.66500 39.68416

# The best period to buy a sunchair is November-February, with an expected price
# of 38.665 euros