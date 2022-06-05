######
#PERMUTATION
######

##### Mean among groups ####
rm(list = ls())
purity <- read.table('purity.txt')
head(purity)

est <- function (i, purity) {
  set.seed(123)
  B <- 100000 # Number of permutations
  T_stat <- numeric(B) # Vector where we will store the values of T*
  n1 <- 8
  n <- 16
  t1 <- purity[1:8,i]
  t2 <- purity[9:16,i]
  T0 <- mean(t2)-mean(t1)
  print(T0)
  for(perm in 1:B){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- purity[permutation,i]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- mean(x2_perm) - mean(x1_perm)
  }
  
  # Permutational distribution of T
  hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
  abline(v=T0,col=3,lwd=2)
  
  plot(ecdf(T_stat))
  abline(v=T0,col=3,lwd=2)
  
  # p-value
  p_val <- sum(T_stat>=T0)/B
  p_val
  return (p_val)
}

p <- c()
for (i in 1:10) {
  p <- cbind(p, est(i, purity = purity))
}
#imposing that the expected proportion  of mutations erroneously judged to be influential 
#among all those judged to be influential is lower than or equal to 1%
p.adjust(p, 'fdr')   #BH
#imposing a probability of at most 1% that at least one of the non-influential mutations is judged as influential
p.adjust(p, 'bonferroni')

### IF ONE GROUP
set.seed(321)
it <- 5000
perm_test <- function(f,n,it=5000){
  i = 1
  stat = rep(NA,it)
  for (i in 1:it){
    n1 <- sample(1:n-1,1)
    g <- f
    perm <- sample(1:length(g))
    g <- g[perm]
    g[1:n1] <- -g[1:n1]
    g[n1+1:length(g)] <- g[(n1+1):length(g)]
    stat[i] <- mean(g)
  }
  return(stat)
}

p_val <- rep(NA,length(df[1,]))
true_stat <- rep(NA,length(df[1,]))

for (k in 1:(length(df[1,]))){
  true_stat[k] <- mean(df[,k])
  stat <- perm_test(df[,k],length(df[,k]),it)
  
  par(mfrow=c(1,2))
  hist(stat,breaks=10)
  abline(v=true_stat[k],col=3,lwd=2)
  plot(ecdf(stat))
  abline(v=true_stat[k],col=3,lwd=2)
  
  p_val[k] <- sum(stat>=true_stat[k])/it
}
true_stat
p_val
cbind(true_stat, p_val)
