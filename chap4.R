# chap4

# 4.1 Kruskal-Wallis 秩和检验
# 这里直接用大样本近似吧...至于精确解，书上的代码嵌套了至少5层循环，
# 感觉可读性不大了...所以就暂时不写了

# 打样本近似[卡方]
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/wtloss.txt")

kruskal.test(d[, 1], d[, 2])


# 4.2 正态计分检验
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/wtloss.txt")
d <- d[order(d[, 1]), ]
n1 <- sum(d[, 2] == 1)
n2 <- sum(d[, 2] == 2)
n3 <- sum(d[, 2] == 3)
n <- nrow(d)
r <- rank(d[, 1], ties.method = "first")  # 按照书上表格这里应为first
w <- qnorm(r/(n+1))
z <- cbind(d, r, w)  # z 即为上面的表格

# T, p值的计算
nn <- sum(sum(w[z[, 2] == 1])^2/n1, sum(w[z[,2]==2])^2/n2, sum(w[z[, 2]==3])^2/n3)
T <- (n-1)*nn/sum(w^2)
pchisq(T, 3-1, low=FALSE)

# 4.3 Jonckheere-Terpstra检验
library(clinfun)
jonckheere.test(d[, 1], d[, 2], alternative = "increasing")


# 4.5 Friedman检验
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/blead.TXT")
friedman.test(as.matrix(d))

# 4.6 Kendall协同系数检验
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/blead.TXT")
R <- apply(d, 2, sum)
m <- nrow(d)
n <- ncol(d)
S <- sum((R-m*(n+1)/2)^2)
W <- 12*S/m^2/(n^3-n)
pchisq(m*(n-1)*W, n-1, lower.tail = F)

# 4.7 完全区组设计：Cochran检验

# way1:大样本近似
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/candid.TXT")
n <- apply(d, 2, sum)
N <- sum(n)
L <- apply(d, 1, sum)
k <- dim(d)[2]
Q <- (k*(k-1)*sum((n-mean(n))^2))/(k*N-sum(L^2))
pvalue <- pchisq(Q, k-1, lower.tail = F)

# way2: RVAideMemoire
library(RVAideMemoire)
response <- unlist(data.frame(t(d)))
fact <- gl(4,1,80,labels=LETTERS[1:4])
block <- gl(20,4,labels=letters[1:20])
cochran.qtest(response~fact|block)

# 4.8 完全区组设计：page检验
# Way1: concord
library(concord)
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/blead1.TXT")
page.trend.test(d)

# Way2:大样本近似
rd <- apply(d, 1, rank)
R <- apply(rd, 1, sum)
L <- sum(R*1:length(R))
k <- dim(d)[2]
b <- dim(d)[1]
m <- b*k*(k+1)^2/4
s <- sqrt(b*(k^3-k)^2/144/(k-1))
Z <- (L-m)/s
pvalue <- pnorm(Z, lower.tail = F)


# 4.9 不完全区组设计：Durbin检验

# 没有打结的情况
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/mater.TXT")
k <- max(d[, 2])
b <- max(d[, 3])
t <- length(d[d[,3]==1,1])
r <- length(d[d[,2]==1,1])
R <- d
for( i in 1:b) R[d[,3]==i,1]=rank(d[d[,3]==i,1])
RV <- NULL
for(i in 1:k) RV = c(RV, sum(R[R[,2]==i,1]))
D <- 12*(k-1)/(r*k*(t^2-1))*sum((RV-r*(t+1)/2)^2)
p.value <- pchisq(D, k-1, lower.tail = F)

# 打结的情况
A <- sum(R[,1]^2)
C <- b*t*(t+1)^2/4
D <- (k-1)*sum((RV-r*(t+1)/2)^2)/(A-C)
p.value <- pchisq(D, k-1, lower.tail = F)


# 课后习题
# T1
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/4.10.1.TXT")

# K-W
kruskal.test(d[, 1], d[, 2])

# J-T
# 按照中位数从小到大的顺序重新排列数据
medians <- aggregate(V1~V2, d, median)
orders <- order(medians$V1)  # 5 1 2 4 3
rep.orders <- NULL
for(i in orders){rep.orders <- c(rep.orders, 
                                 rep(i, nrow(d[d$V2==i, ])))}
ordered.d <- d[match(rep.orders, d$V2), ]

library(clinfun)
# H1: theta5 <= theta1 <= ... <= theta3
jonckheere.test(ordered.d[, 1], ordered.d[, 2], alternative = "increasing")


#T2
d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/4.10.2.txt")
library(clinfun)
jonckheere.test(d[, 1], d[, 2], alternative = "increasing")

#T9
# d <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/4.10.9.txt")
d <- read.table(file.choose())                                                                  
k <- max(d[, 2])
b <- max(d[, 3])        
t <- length(d[d[,3]==1,1])
r <- length(d[d[,2]==1,1])
R <- d
for( i in 1:b) R[d[,3]==i,1]=rank(d[d[,3]==i,1])
RV <- NULL
for(i in 1:k) RV = c(RV, sum(R[R[,2]==i,1]))
D <- 12*(k-1)/(r*k*(t^2-1))*sum((RV-r*(t+1)/2)^2)
p.value <- pchisq(D, k-1, lower.tail = F)
