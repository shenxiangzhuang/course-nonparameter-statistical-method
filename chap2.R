# 章节内容

# 2.1.1 获取符号检验结果
sign.test <- function(x, p, q0){
  s1 = sum(x < q0)
  s2 = sum(x > q0)
  n = s1 + s2
  p1 = pbinom(s1, n, p)
  #print(p1)
  p2 = 1 - pbinom(s1-1, n, p)
  #print(p2)
  # 书上这里可能有错误，下面符号搞反了，这里已更正
  if(p1 > p2) m1 = "One tail test: H1: Q < q0"
  else m1 = "One tail test: H1: Q > q0"
  p.value = min(p1, p2)
  m2 = "Two tails test"
  p.value2 = 2 * p.value
  if(q0 == median(x)){p.value=0.5; p.value2=1}
  list(Sign.test1=m1, p.values.of.one.tail.test=p.value, 
       p.value.of.two.tail.test=p.value2)
}


library(readr)
ExpensiveCities <- read_csv("~/Documents/MyPrograming/R/非参数统计/data/ExpensiveCities.csv")
#View(ExpensiveCities)
hist(ExpensiveCities$`Index_including rent`)
library(DescTools)
SignTest(ExpensiveCities$`Index_including rent`, mu=64, alternative="greater")
sign.test(ExpensiveCities$`Index_including rent`, 0.5, 64)


# 2.1.2 
# 基于符号检验的中位数及分位点的置信区间

# way1
mci <-  function(x, alpha=0.05){
  x = sort(x)  # 要先排好序
  n = length(x)
  b = 0
  i = 0
  while((b<=alpha/2) & (i<=floor(n/2))){
    b = pbinom(i, n, .5)
    k1 = i  # 从前往后
    k2 = n-i+1  #从后往前
    a = 2*pbinom(k1-1, n, .5)
    i = i+1;
  }
  
  z = c(k1, k2, a, 1-a)
  #print(z)
  z2 = "Entire range!"
  if(k1>1)
    out = list(Confidence.level=1-a, CI=c(x[k1], x[k2]))
  else
    out = list(Confidence.level=1-2*pbinom(0, n, 0.5), CI=z2)
  out
}


# way1 test
library(readr)
data2.2 <- read_table("~/Documents/MyPrograming/R/非参数统计/data/tax.txt",
                      col_names = FALSE)
# 注意前面读入的数据格式，这里要索引出来纯数据
mci(data2.2, 1-0.999989)


# way2
mci2 <- function(x, alpha=0){
  x = sort(x)
  n = length(x)
  q = .5
  m = floor(n*q)
  s1 = pbinom(0:m, n, q)
  s2 = pbinom(m:(n-1), n, q, lower.tail = F)
  ss = c(s1, s2)
  nn = length(ss)
  a = NULL
  for(i in 0:m)
  {
    b1 = ss[i+1]
    b2 = ss[nn-i]
    b = b1 + b2
    d = 1 - b
    if(b>1) break
    a = rbind(a, c(b, d, x[i+1], x[n-i]))
  }
  # 书上代码没有体现alpha为缺省值时的输出情况，这里补上
  if (alpha == 0){out = a}  
  else if(a[1, 1]>alpha) {out="alpha is too small, CI=All range"}
  else{
    for(i in 1:nrow(a))
      if(a[i,1]>alpha){out = a[i-1, ];break}
  }
  out
}

# way2 test
# mci2(data2.2)
mci2(data2.2, alpha = 1-0.999989)



# 分位点Q的置信区间

qci <- function(x, alpha=0.05, q=0.25){
  x <- sort(x)
  n = length(x)
  a = alpha/2
  r = qbinom(a, n, q);
  s = qbinom(1-a, n, q);
  CI = pbinom(s, n, q) - pbinom(r-1, n, q)
  if(r == 0) lo <- NA else lo <- x[r]
  if(s == n) up <- NA else up <- x[s+1]
  list(c("Lower limit"=lo, "Upper limit"=up, "1-alpha"=1-alpha,
         "True conf"=CI))
}

# qci test
qci(data2.2, alpha = 1-0.999989, q = 0.5)


# 2.2 Wilcoxon符号秩检验

library(haven)
euroalc <- read_sav("~/Documents/MyPrograming/R/非参数统计/data/euroalc.sav")
euroalc <- sort(euroalc$y)

# 右侧秩和检验
wilcox.test(euroalc-8, alternative = "greater")
# 左侧秩和检验
wilcox.test(euroalc-12.5, alternative = "less")

# 基于Wilcoxon秩和检验的点估计与区间估计[Walsh平均]

WalshConf <- function(x, alpha=0.05){
  x <- sort(x)
  n <- length(x)
  N <- n*(n+1)/2
  
  walsh <- NULL
  for (i in 1:n)
    for(j in i:n)
    {
      walsh <- c(walsh, (x[i]+x[j])/2)
    }
  walsh <- sort(walsh)
  k <- qsignrank(alpha/2, n)
  CI <- c(walsh[k+1], walsh[N-k])
  out <- list(c("中位数的置信区间："=CI, "显著性水平"=alpha))
  out
}

WalshConf(euroalc, 0.05)



# 2.3 正态计分检验[书上有错]
ns <- function(x, m0){
  x1 <- x-m0
  n <- length(x)
  r = rank(abs(x1))
  s = qnorm(.5*(1+r/(n+1)))*sign(x1)
  tt = sum(s)/sqrt(sum(s^2))
  list(pvalue.2sided=2*min(pnorm(tt), pnorm(tt, low=F)), 
       T = tt, s=s, Sn=sum(s))
}

ns(euroalc, 12.5)


# 2.4 Cox-Stuart趋势检验

#参考[https://www.r-bloggers.com/trend-analysis-with-the-cox-stuart-test-in-r/]

library(readr)
TJAir <- read_table("~/Documents/MyPrograming/R/非参数统计/data/TJAir.txt",
                      col_names = FALSE)

TJAir <- as.vector(t(TJAir))  # 按行合并成vector[按列的话就用unlist]

cox.stuart.test <- function(x)
  {
    method = "Cox-Stuart test for trend analysis"
    leng = length(x)
    apross = round(leng) %% 2
    if (apross == 1) {
      delete = (length(x)+1)/2
      x = x[ -delete ] 
    }
    half = length(x)/2
    x1 = x[1:half]
    x2 = x[(half+1):(length(x))]
    difference = x1-x2
    signs = sign(difference)
    signcorr = signs[signs != 0]
    pos = signs[signs>0]
    neg = signs[signs<0]
    
    if (length(pos) < length(neg)) {
      prop = pbinom(length(pos), length(signcorr), 0.5)
      names(prop) = "Increasing trend, p-value"
      rval <- list(method = method, statistic = prop)
      class(rval) = "htest"
      return(rval)
    }
    else {
      prop = pbinom(length(neg), length(signcorr), 0.5)
      names(prop) = "Decreasing trend, p-value"
      rval <- list(method = method, statistic = prop)
      class(rval) = "htest"
      return(rval)
    }
  }

cox.stuart.test(TJAir)


# 2.5 关于随机性的游程检验

# way1
# 注意tseries的run.test仅能用于近似正态的游程检验
library(tseries)
library(haven)
run02 <- read_sav("~/Documents/MyPrograming/R/非参数统计/data/run02.sav")
run02 <- as.vector(t(run02))
y <- factor(sign(run02-median(run02)))
tseries::runs.test(y)



# way2[书上的实现]

run.test <- function(y, cut=0){
  
  # 如果不是二元数据，指定参数cut为中位数可以进行转化
  if(cut!=0)
    x = (y > cut)*1
  else
    x = y
  
  N = length(x)
  
  # 计算游程数
  k = 1
  for(i in 1:(N-1))
    if (x[i] != x[i+1]) k = k + 1
  r = k
  m = sum(1-x)
  n = N-m
  
  P1 = function(m, n, k){
    2*choose(m-1, k-1)*choose(n-1, k-1)/choose(m+n, n)
  }
  
  P2 = function(m, n, k){
    choose(m-1, k-1)*choose(n-1, k)/choose(m+n, n)
    + choose(m-1, k)*choose(n-1, k-1)/choose(m+n, n)
  }
  
  r2 = floor(r/2)
  if(r2 == r/2){
    pv = 0
    for(i in 1:r2)
      pv = pv + P1(m, n, i)
    for(i in 1:(r2-1))
      pv = pv + P2(m, n, i)
  }
  else{
    pv = 0
    for(i in 1:r2)
      pv = pv + P1(m, n, i)
    for(i in 1:r2)
      pv = pv + P2(m, n, i)
  }
  
  if(r2 == r/2)
    pv1 = 1 - pv + P1(m, n, r2)
  else
    pv1 = 1- pv + P2(m, n, r2)
  
  z = (r-2*m*n/N-1)/sqrt(2*m*n*(2*m*n-m-n)/(m+n)^2/(m+n-1))
  ap1 = pnorm(z)
  ap2 = 1 - ap1
  tpv = min(pv, pv1)*2
  list(m=m, n=n, R=r, 
       Exact.pvalue1=pv, Exact.pvalue2=pv1, Eaxct.2sided.pvalue=tpv,
       Aprox.pvalue1 = ap1, Aprox.pvalue2=ap2, Aprox.2sided.pvalue=min(ap1, ap2)*2)
        
}

run.test(run02, median(run02))



# 习题2.6

# T1

# S1: Load data
library(readr)
data2.6.1 <- read.delim("~/Documents/MyPrograming/R/非参数统计/data/2.6.1.TXT", header = FALSE)
data2.6.1 <- as.vector(t(data2.6.1))
# S2: 计算
# median
data.median <- median(data2.6.1)
# qci 

qci <- function(x, alpha=0.05, q=0.25){
  x <- sort(x)
  n = length(x)
  a = alpha/2
  r = qbinom(a, n, q);
  s = qbinom(1-a, n, q);
  CI = pbinom(s, n, q) - pbinom(r-1, n, q)
  if(r == 0) lo <- NA else lo <- x[r]
  if(s == n) up <- NA else up <- x[s+1]
  list(c("Lower limit"=lo, "Upper limit"=up, "1-alpha"=1-alpha,
         "True conf"=CI))
}

# q = 0.25
qci(data2.6.1, alpha = 0.05, q = 0.25)
# q = 0.75
qci(data2.6.1, alpha = 0.05, q = 0.75)

# H0:median=1200
sign.test <- function(x, p, q0){
  s1 = sum(x < q0)
  s2 = sum(x > q0)
  n = s1 + s2
  p1 = pbinom(s1, n, p)
  #print(p1)
  p2 = 1 - pbinom(s1-1, n, p)
  #print(p2)
  # 书上这里可能有错误，下面符号搞反了，这里已更正
  if(p1 > p2) m1 = "One tail test: H1: Q < q0"
  else m1 = "One tail test: H1: Q > q0"
  p.value = min(p1, p2)
  m2 = "Two tails test"
  p.value2 = 2 * p.value
  if(q0 == median(x)){p.value=0.5; p.value2=1}
  list(Sign.test1=m1, p.values.of.one.tail.test=p.value, 
       p.value.of.two.tail.test=p.value2)
}

sign.test(data2.6.1, 0.5, 1200)  # p=0.0153<0.05,拒绝H0，中间值不在1200
sign.test(data2.6.1, 0.25, 750)  # p=0.16>0.05,没有理由拒绝H0，认为下四分位不少于750


# T7
library(readr)
data2.6.7 <- read.delim("~/Documents/MyPrograming/R/非参数统计/data/2.6.7.TXT", header = FALSE)
data2.6.7 <- as.vector(t(data2.6.7))
tseries::runs.test(factor(data2.6.7))  # p=0.7123>0.05,所以接受H0，认为其随机


