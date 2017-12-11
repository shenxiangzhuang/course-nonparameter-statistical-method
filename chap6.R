# chap6

# 6.1 Spearman
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/DM1.txt")
x <- d[, 2]
y <- d[, 1]
rx <- rank(x)
ry <- rank(y)
rsd <- rbind(rx, ry, (rx-ry)^2)
cor.test(x, y, method = "spearman", exact = T)

# 6.2 Kendall 相关检验
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/CPIESI.TXT")
x <- d[, 2]
y <- d[, 1]
cor.test(x, y, method = "kendall", exact = T)

# 6.3 Goodman-Kruskal's gamma

# way1: only get gamma and confident interval
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/incsat.txt")
M <- matrix(nrow = length(unique(d$V1)), ncol = length(unique(d$V2)))
index <- 1
for(i in 1:length(unique(d$V1))){
  for(j in 1:length(unique(d$V2))){
    M[i, j] <- d$V3[index]
    index <-  index + 1
    #print(index)
  }
}
library(DescTools)
GoodmanKruskalGamma(M, direction="column", conf.level=0.95)

# way2: get gamma, pvalue and ASE
xx <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/incsat.txt")
x = xx[, 1]
y = xx[, 2]
w = xx[, 3]
n1 = max(x)
n2 = max(y)
WW <- matrix(w, byrow=T, nrow=n1)
DD = CC = matrix(0, nrow = n1, ncol = n2)
for (i in 1:n1) {
  for(j in 1:n2){
    CC[i, j] = sum((x>i)*(y>j)*w) + sum((x<i)*(y<j)*w)
    DD[i, j] = sum((x>i)*(y<j)*w) + sum((x<i)*(y>j)*w)
  }
}

nc <- sum(WW*CC)/2
nd <- sum(WW*DD)/2
G <- (nc-nd)/(nc+nd)
ASE <- 1/(nc+nd)^2*sqrt(sum(WW*(2*nd*CC-2*nc*DD)^2))
pvalue <- 2*(1-pnorm(G/ASE))
CI95 <- c(G-1.96*ASE, G+1.96*ASE)
list(G=G, ASE=ASE, CI95=CI95, pvalue=pvalue)


# 6.4 Somers'd 相关检验
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/incsat.txt")
M <- matrix(nrow = length(unique(d$V1)), ncol = length(unique(d$V2)))
index <- 1
for(i in 1:length(unique(d$V1))){
  for(j in 1:length(unique(d$V2))){
    M[i, j] <- d$V3[index]
    index <-  index + 1
    #print(index)
  }
}
library(DescTools)
# Somers' D C|R    (or Y|X)
SomersDelta(M, direction="column", conf.level=0.95)
# Somers' D R|C    (or X|Y)
SomersDelta(M, direction="row", conf.level=0.95)


# Theil非参数回归和几种稳健回归
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/CPIGINI.txt",
                header = T)

library(mblm)
# Theil
fit1 <- mblm(GINI~CPI, d, repeated = F)
# Siegel
fit2 <- mblm(GINI~CPI, d, repeated = T)

summary(fit2)
anova(fit2)
confint(fit2)

# LMS. LTS, S
library(MASS)
Y = d$GINI
X = d$CPI

lms = lqs(Y~X, method = "lms")
lts = lqs(Y~X, method = "lts")
S = lqs(Y~X, method = "S")
summary(S)

