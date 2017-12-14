# Chap7

# 7.1 Kolmogorov-Smirnov单样本检验及一些正态检验
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/ind.TXT")
# Marsaglia 2003
ks.test(x, "pnorm", 15, 0.2, exact = T)  # p-value = 0.0147

library(nortest)
# Lilliefors 1967 && 还可以进行其他很多检验，这里不列举了
lillie.test(x$V1)  # p-value = 0.6847

# 7.2 Kolmogorov-Smirnov两样本检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/ks2.txt",
                header = F)
x <- z[z[, 2]==1, 1]
y <- z[z[, 2]==2, 1]
# way1
ks.test(x, y)
# way2
library(fBasics)
ks2Test(x, y)

# 7.3 Pearson Chi-Square

# Poisson distribution test
Ob <- c(490, 334, 68, 16)
n <- sum(Ob)
lambda <- t(0:3)%*%Ob/n  # %*% Matrix Multiplication
p <- exp(-lambda)*lambda^(0:3)/factorial(0:3)
E <- p*n
Q <- sum((E-Ob)^2/E)
pvalue = pchisq(Q, 2, low=F)

# 零假设为多项式分布的时候
chisq.test(Ob, p=p, rescale.p = T)
