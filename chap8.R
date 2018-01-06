# 列联表

# 8.1 二维列联表的齐性检验和独立性的卡方检验
# 齐行检验
y <- matrix(scan("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/wid.TXT"), 3, 2, b=T)
chisq.test(y)  # p-value = 0.5839

# 独立性检验
y <- matrix(scan("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/shop.TXT"), 3, 3, b=T)
chisq.test(y)  # p-value = 0.0009203

# 8.2 低纬列联表的Fisher精确检验
# way1:针对表格形式的数据[频数方阵]
x = read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/stroke.txt", header = F)
fisher.test(x)  # p-value = 0.002242

# way2:针对不同组合的频数数据
x = read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/strokeA.txt", header = T)
M = xtabs(X35~., x)  # 这里书上用的Freq，其实就是指频数列，这里为X35
fisher.test(M)  # p-value = 0.001759


# 8.3 两个比例的比较[相对风险和胜算比]
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/stroke.txt")
p1 <- x[1, 1]/sum(x[1, ])
p2 <- x[2, 1]/sum(x[2, ])
# 两比例之差
pdif1 <- p1 - p2
se1 <- sqrt(p1*(1-p1)/sum(x[1, ])+p2*(1-p2)/sum(x[2, ]))
pdifc1 <- c(p1-p2-1.96*se1, p1-p2+1.96*se1)
# 相对风险
rr1 <- p1/p2
ser1 <- sqrt((1-p1)/x[1,1]+(1-p2)/x[2,1])
rrc1 <- c(rr1*exp(-1.96*ser1), rr1*exp(1.96*ser1))
# 胜算比
or1 <- (p1/(1-p1))/(p2/(1-p2))
seor1 <- sqrt(sum(1/x))
orc1 <- c(or1*exp(-1.96*seor1), or1*exp(1.96*seor1))
list(dif=pdif1, difCI=pdifc1, RR=rr1, RRCI=rrc1, ORCI=orc1)

# 8.4 Cochran-Mantel-Haenszel估计
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/hospital.txt")
tmp <- array(c(x[, 4]), dim=c(2,2,4), dimnames=list(effect=c("Y", "N"),
                                                   med=c("A", "B"), 
                                                   hospital=c("I", "II", "III", "IV")))

tab <- ftable(.~med+effect, tmp)  # 可以看作依med和effect作为分类的“自变量”
mantelhaen.test(tmp)


# 8.6对数线性模型与高纬列联表的独立性检验简介
x <- read.table("/home/shensir/Documents/MyPrograming/R/NonparameterStat-Course-WIth-R/data/wmq.TXT", header = T)
xt <- xtabs(Count~., x)
a <- loglin(xt, list(1:2, c(1:3)), param = T)
