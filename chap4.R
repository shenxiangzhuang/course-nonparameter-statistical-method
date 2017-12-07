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

# 4.8

