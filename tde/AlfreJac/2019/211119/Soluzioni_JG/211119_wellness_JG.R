wellness <- read.table("wellness.txt", header = T)
n <- dim(wellness)[1]
p <- dim(wellness)[2]
p_val <- rep(NULL, p)
T_0 <- rep(NULL, p)
B <- 5000
#Individual test:
#H0: mean_i < 0 vs H1
#Test statistic: T_0 = sample_mean_i, reject H0 when it is high

for(i in 1:p){
  x <- wellness[,i]
  T_0[i] = mean(x)
  set.seed(321)
  T_stat <- rep(NULL, B)
  for(j in 1:B){
    #Idea of the permutation:
    #for each patient we try and attribute the wellness result of another random treatment
    #So for instance a possible permutation could be:
    # patient 1 -> his/her result of treatement 5
    # patient 2 -> his/her result of treatement 9
    #...
    #So we compare the overall result of the treatement with the possible result of the effects
    #of random treatement on each individual patient
    permutation <- sample(1:p, size = n, replace = TRUE)
    x_perm <- rep(NULL, n)
    for(k in 1:n){
      x_perm[k] <- wellness[k, permutation[k]]
    }
    T_stat[j] = mean(x_perm)
  }
  p_val[i] <- sum(T_stat >=T_0[i])/B
}
p_val
#So the treatements with a meaningful impact are treatement 4 and 10
#makes sense, considering:
colMeans(wellness)
which(colMeans(wellness)>0.1)

#we want to control the FWER, Bonferroni does just that:
alphaB <- 0.1/11
which(p_val<alphaB)
