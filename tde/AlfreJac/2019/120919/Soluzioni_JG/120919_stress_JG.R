stress<-read.table("stress.txt", header = T)
n <- dim(stress)[1]
p <- dim(stress)[2]
mod <- c(rep(1, 10), rep(0,5))
modif <- stress[which(mod==1),]
unmodif <- stress[which(mod==0),]
B <- 10000
nm <- dim(modif)[1]
nu <- dim(unmodif)[1]

T_0 <- rep(NULL, p)
p_val <- rep(NULL, p)
#Test:
#H0: median_mod <= median_unmod vs H1
#T_stat median_mod - median_unmod -> reject when high
#we permute over the realization of the results of the stress test for each unit, regardless of wether it was modified or not

for(i in 1:p){
  T_0[i] <- median(modif[,i]) - median(unmodif[,i])
  x_base <- stress[,i]
  set.seed(123)
  T_stat <- rep(NULL, B)
  for(j in 1:B){
    perm <- sample(1:n)
    x_perm <- x_base[perm]
    T_stat[j] <- median(x_perm[1:nm]) - median(x_perm[(nm+1):n])
  }
  p_val[i] <- sum(T_stat >= T_0[i])/B
}
p_val
T_0

#We need to control FDR, using Benjamini &Kekuteli:
alpha <- 0.25
P_ord <- p_val[order(p_val, decreasing = FALSE)]
z <- 1:length(P_ord)
zz <- P_ord<=z*alpha/length(P_ord)
zz
m <- 4
rej <- order(p_val, decreasing = FALSE)[1:m]
rej #we reject H0 for those stress tests -> in here we have an improvement
p_val[rej]
