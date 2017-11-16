# 第三章

# 3.1 两样本和多样本的Brown-Mood中位数检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/
                NonparameterStat-Course-WIth-R/data/
                salary.TXT")  # 第二列为对应分组

k <- unique(z[, 2])  # 组别
m <- median(z[, 1])  # 两组混合的，工资的，中位数

m1 <- NULL;  # a b
m2 <- NULL;  # m-a n-b

for (i in k){
  m1 <- c(m1, sum(z[z[, 2] == i, 1] > m))  # a b
  m2 <- c(m2, sum(z[z[, 2] == i, 1] <= m))  # m-a n-b
}

C <- rbind(m1, m2)  # 2x2矩阵

# 我们有三种办法来计算Brown-Mood检验的P值

# way1：Ｃ为Fisher精确检验的特殊形式
pVal <- fisher.test(C, alternative = "less")  # 左侧检验 (0.07781)

# way2: 精确地直接使用超几何分布
# 为了方便，先提取出a, b, m, n
a <- C[1, 1]
b <- C[1, 2]
m <- C[2, 1]+a
n <- C[2, 2]+b
pValLeft <- phyper(a, m, n, a+b)  # left (0.07780674)
pValRight <- phyper(b, n, m, a+b)  # Right

# way3: 大样本近似
t <- a+b
# 同样以左侧检验测试，其中0.5是用于连续性修正
pValAsy <- pnorm((a+0.5-m*t/(m+n))/sqrt(m*n*t*(m+n-t)/(m+n)^3))  # (0.07824383)


# 3.2 Wilcoxon(Mann-Whitney)秩和检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/
                NonparameterStat-Course-WIth-R/data/
                salary.TXT")  # 第二列为对应分组

# 这里也是两种不同的思路

# way1: wilcox.test(注意此时的p值基于大样本的正态近似来计算的)
df <- data.frame(z)
x <- df[1][df[2] == 1]
y <- df[1][df[2] == 2]
pVal <- wilcox.test(x, y, alternative = "less")  # 0.01352

# way2: 基于Wxy(Wyx)
# 依旧以左侧检验为例，选取Wyx
m <- length(x)
n <- length(y)
Wyx <- sum(outer(x, y, "-")>0)
pVal <- pwilcox(Wyx, n, m)  # (0.01352166)

# Mx-My的点估计与区间估计

pointEst <- median(outer(x, y, "-"))
D <- sort(as.vector(outer(x,y, "-")))
Alpha <- 0.05
wHalfAlpha <- qwilcox(Alpha/2, m, n) 

lowIndex <- wHalfAlpha
upIndex <- m*n+1-wHalfAlpha
confInv <- c(D[lowIndex], D[upIndex])  # -3916  -263


# 3.3 正态计分检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/
                NonparameterStat-Course-WIth-R/data/
                salary.TXT")  # 第二列为对应分组

df <- data.frame(z)
x <- df[1][df[2] == 1]
y <- df[1][df[2] == 2]

w <- cbind(c(x ,y), c(rep(1,17), rep(2, 15)))  # 教材这句少了一个括号
w <- w[order(w[, 1]), ]
w <- cbind(w, 1:32, qnorm((1:32)/(17+15+1)))

T <- sum(w[w[,2] == 1, 4])
w2 <- sum(w[, 4]^2)
S <- sqrt(m*n*w2/(m+n-1)/(m+n))
Z <- T/S

pVal <- pnorm(Z)  # 0.01925286



