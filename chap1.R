# 章节内容

# 图1.1
x <- rnorm(500)
y <- rexp(500, 1)
par(mfrow=c(2, 3))
hist(x, mian="a. Histogram of x")
boxplot(x, main="b. Boxplot of x")
qqnorm(x, main="c. Normal Q-Q Plot of x")

hist(y, mian="a. Histogram of y")
boxplot(y, main="b. Boxplot of y")
qqnorm(y, main="c. Normal Q-Q Plot of y")



# 图1.2
# 书中写的可能有误[不是标准正态，mean=1]
BoxCox <- function(x, lambda){
  if (lambda == 0){
    return(log(x))
  }
  return((x**lambda-1)/lambda)
}


x <- exp(rnorm(5000, 0, 1))
lambdas <- c(1.5, 1, 0.75, 0.4, 0, -0.25, -0.5, -0.75)
par(mfrow=c(2, 4))
for (lambda_index in 1:length(lambdas)){
  y <- BoxCox(x, lambdas[lambda_index])
  hist(y, main = substitute(paste(lambda, " = ", lambda_value),
                            list(lambda_value=lambdas[lambda_index])))
}



# 课后题

# T3(1)
x1 <- rnorm(100, 0, 1)
x2 <- rnorm(20, 3, 3)
x12 <- c(x1, x2)
par(mfrow(1, 3))
hist(x12)
boxplot(x12)
qqnorm(x12)

# T3(2)
t.test(x12, mu=0, alternative = "greater")
# t = 3.3857, df = 119, p-value = 0.0004813, 显著，拒绝原假设，均值不为零[大于0]

# T3(3)
t.test(x1, mu=0, alternative = "greater")
# t = 1.0365, df = 99, p-value = 0.1513, 不显著，没有理由拒绝原假设


# T4

data <- rnorm(100, 20, 1)
lambdas <- c(1.5, 1, 0.75, 0.4, 0, -0.25, -0.5, -0.75)
par(mfrow=c(2, 4))
for (lambda_index in 1:length(lambdas)){
  y <- BoxCox(data, lambdas[lambda_index])
  hist(y, main = substitute(paste(lambda, " = ", lambda_value),
                            list(lambda_value=lambdas[lambda_index])))
}



