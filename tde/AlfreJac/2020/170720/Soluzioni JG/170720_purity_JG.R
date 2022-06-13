purity <- read.table("purity.txt", header = T)
old <- purity[1:8, ]
new <- purity[9:16, ]
p <- dim(purity)[2]
alpha <- 0.1

n1 <- dim(new)[1]
n2 <- dim(old)[1]
n <- n1 + n2
p_val <- rep(NULL, p)


for(i in 1:p){
  x_pooled <- c(new[,i],old[,i])

  T0[i] <- mean(new[,i]) - mean(old[,i])     #we want to reject H0: old is better than new, i.e reject when T0 is large
  T0[i]
  set.seed(123)
  B <- 5000 # Number of permutations
  T_stat <- numeric(B) # Vector where we will store the values of T*
  
  for(perm in 1:B){
    # permutation:
    permutation <- sample(1:n)
    x_perm <- x_pooled[permutation]
    x1_perm <- x_perm[1:n1]
    x2_perm <- x_perm[(n1+1):n]
    # test statistic:
    T_stat[perm] <- (mean(x1_perm) - mean(x2_perm))
  }
  
  p_val[i] <- sum(T_stat>=T0[i])/B
}

T0
p_val

#we can assume independence on the p_values, being defined by looking at different purity parameters
#and computing permutation on those
#So we control FDR with Benjamini & Yekuteli
P_ord <- p_val[order(p_val, decreasing = FALSE)]
z <- 1:length(P_ord)
zz <- P_ord<=z*alpha/length(P_ord)
zz
m <- 5
rej <- order(p_val, decreasing = FALSE)[1:m]
rej #we reject H0 for those purity parameters
p_val[rej]

#So for the purity parameters indexed by rej the new dealer is significantly better than the old one

#request c is a request of level 1% over the Family Wise Error Rate, controlled by Bonferroni
alphaB <- 0.01/10
which(p_val < alphaB)
