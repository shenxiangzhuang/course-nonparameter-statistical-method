# 第三章

# 3.1 两样本和多样本的Brown-Mood中位数检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/salary.TXT")  # 第二列为对应分组

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

# way1: wilcox.test(注意此时的p值默认基于大样本的正态近似来计算的)
df <- data.frame(z)
x <- df[1][df[2] == 1]
y <- df[1][df[2] == 2]
pVal <- wilcox.test(x, y, alternative = "less")  # 0.01352
#wilcox.test(x, y, exact = T, alternative = "less", conf.int = T, conf.level = 0.95)

# way2: 基于Wxy(Wyx)
# 依旧以左侧检验为例，选取Wyx
m <- length(x)
n <- length(y)
Wyx <- sum(outer(x, y, "-")>0)
pVal <- pwilcox(Wyx, n, m)  # (0.01352166)

# Mx-My的点估计与区间估计[这里仅仅为双边检验的区间估计]

pointEst <- median(outer(x, y, "-"))
D <- sort(as.vector(outer(x,y, "-")))
Alpha <- 0.05
wHalfAlpha <- qwilcox(Alpha/2, m, n) 

lowIndex <- wHalfAlpha
upIndex <- m*n+1-wHalfAlpha
confInv <- c(D[lowIndex], D[upIndex])  # -3916  -263


# 3.3 正态计分检验
z <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/salary.TXT")  # 第二列为对应分组

df <- data.frame(z)
x <- df[1][df[2] == 1]
y <- df[1][df[2] == 2]
m <- length(x)
n <- length(y)

w <- cbind(c(x ,y), c(rep(1,17), rep(2, 15)))  # 教材这句少了一个括号
w <- w[order(w[, 1]), ]
w <- cbind(w, 1:32, qnorm((1:32)/(17+15+1)))

T <- sum(w[w[,2] == 1, 4])
w2 <- sum(w[, 4]^2)
S <- sqrt(m*n*w2/(m+n-1)/(m+n))
Z <- T/S

pVal <- pnorm(Z)  # 0.01925286


# 3.4 成对数据的检验

# 符号检验？？？
# 书上下面用的psignrank， 感觉应该用pbinom，先存疑...
# two-side
pVal <- 2*pbinom(min(sum(x<y), sum(x>y)), n, 1/2)
# H1: MD > MD0
pVal <- pbinom(sum(x<y), n, 1/2)


# Wilcoxon秩和检验
# Way1: build-in funciton to test
data <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/bp.txt")
x <- data$V1
y <- data$V2
n <- length(x)
wilcox.test(x, y, paired = TRUE, alternative = "greater")  # 0.01367

# Way2: step by step
Di <- x - y
Di_rank <- rank(abs(Di))
W_pos <- sum(Di_rank * ((x-y)>0))
W_neg <- sum(Di_rank * ((x-y)<0))
w = min(W_pos, W_neg)
psignrank(w, n)  # 0.01367188


# 3.5 McNemar检验
# 精确计算，暂略

# 大样本近似
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/athletefootp.txt")
x <- x[, -1] # 去除序号列
# 书上这里的代码好像是有错的，与书上的数据不符合...这里换种方法来求
n11 <- length(x[((x[, 1]==1)&(x[, 2]==1)), 1])
n12 <- length(x[((x[, 1]==1)&(x[, 2]==0)), 1])
n21 <- length(x[((x[, 1]==0)&(x[, 2]==1)), 1])
n22 <- length(x[((x[, 1]==0)&(x[, 2]==0)), 1])

McNemar <- (n12-n21)^2/(n12+n21)
pvalue <- 1 - pchisq(McNemar, df=1)
list(McNemar=McNemar, pvaluetwoside=pvalue)


# stat::mcnemar
Performance <-
  matrix(c(n11, n21, n12, n22),
         nrow = 2,
         dimnames = list("A Drug" = c("Good", "Bad"),
                         "B Drug" = c("Good", "Bad")))
Performance
mcnemar.test(Performance, correct = FALSE)  # 0.001091
# continuity correction 
mcnemar.test(Performance, correct = TRUE)   # 0.0022
mcnemar.exact


# 3.6 Cohen's Kappa系数
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/music.txt")

w <- matrix(x[, 3], byrow = TRUE, ncol = 2)
I <- nrow(w)
n <- sum(w)
w <- w/n
pa <- sum(diag(w))
pe <- sum(apply(w, 1, sum)*apply(w, 2, sum))
kap <- (pa-pe)/(1-pe)
A <- sum(diag(w)*(1-(apply(w, 1, sum)+apply(w, 2, sum))*(1-kap))^2)
tempB <- matrix(rep(apply(w, 1, sum), I)+
                  rep(apply(w, 2, sum), each=I), byrow = TRUE, ncol = I)
diag(tempB) <- 0
B <- (1-kap)^2*sum(w*tempB^2)
CC <- (kap-pe*(1-kap))^2
ASE <- sqrt((A+B-CC)/(1-pe)^2/n)
list(kappa=kap, ASE=ASE, CI=c(kap-1.96*ASE, kap+1.96*ASE))


# 3.7习题
# T1 采用wilcoxon秩和检验
data <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.1.txt")
x <- data$V1[data$V2==1]
y <- data$V1[data$V2==2]

wilcox.test(x, y, exact = T, alternative = "less", conf.int = T, conf.level = 0.95)


# T2
# wilcoxon test
data <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.2.txt")
x <- data$V1[data$V2==1]
y <- data$V1[data$V2==2]
# annot compute exact p-value with ties
# p=0.001041, 拒绝H0  D.conf.inv = [4.30004 13.30000]
wilcox.test(x, y, alternative = "two.side", conf.int = T, conf.level = 0.95)


# t test,拒绝H0  D.conf.inv = [3.916182 13.346929]
# p=0.00174
t.test(x, y, alternative = "two.side", conf.level = 0.95)


#T3
# Brown-Mood
z <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.3.txt")  
# 第二列为对应分组
k <- unique(z[, 2])  # 组别
m <- median(z[, 1])  # 两组混合的，工资的，中位数

m1 <- NULL;  # a b
m2 <- NULL;  # m-a n-b

for (i in k){
  m1 <- c(m1, sum(z[z[, 2] == i, 1] > m))  # a b
  m2 <- c(m2, sum(z[z[, 2] == i, 1] <= m))  # m-a n-b
}

C <- rbind(m1, m2)  # 2x2矩阵
# p-value = 0.1556   [0.0489375 1.4470395]
fisher.test(C, alternative = "two.side", conf.int = T, conf.level = 0.95)


# wilcoxon test
data <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.3.txt")
x <- data$V1[data$V2==1]
y <- data$V1[data$V2==2]

# p-value = 0.0004883 conf.inv=[-7.1 -1.6]
wilcox.test(x, y, paired = T,alternative = "two.side", conf.int = T, conf.level = 0.95)


# t test
# p-value = 0.0005278 conf.inv=[ -6.465888 -2.417446]
t.test(x, y, paired = T, alternative = "two.side", conf.level = 0.95)


# T4
# Brown-Mood
z <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.4.txt")  
# 第1列为对应分组
k <- unique(z[, 1])  # 组别
m <- median(z[, 2])  # 两组混合的，工资的，中位数

m1 <- NULL;  # a b
m2 <- NULL;  # m-a n-b

for (i in k){
  m1 <- c(m1, sum(z[z[, 1] == i, 2] > m))  # a b
  m2 <- c(m2, sum(z[z[, 1] == i, 2] <= m))  # m-a n-b
}

C <- rbind(m1, m2)  # 2x2矩阵
# p-value = 0.9714 can not reject H0
fisher.test(C, alternative = "less")


# wilcoxon 
data <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/3.7.4.txt")
A <- data$V2[data$V1==1]
B <- data$V2[data$V1==2]
wilcox.test(A, B, "less" )  # p-value = 0.9087, can not reject H0


m <- prophet(df[1:(length(df$ds)-9), ])
future <- make_future_dataframe(m, periods = 9, freq = 'm',)
tail(future,9)


