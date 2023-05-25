rm(list=ls())
require(agridat)
require(reshape2)
data(kang.peanut)
set.seed(42)
Y0 <- acast(kang.peanut, gen~env, value.var="yield", fun.aggregate = mean)


# Specify the significance level and the number of bootstrap samples

alpha <- 0.05

N.Boot <- 1000


# Program start (no changes are needed below)

n <- nrow(Y0)
p <- ncol(Y0)
M <- min(n - 1, p)


# Variables are centred, but not scaled

Y <- scale(Y0, center = TRUE, scale = FALSE)


# Singular value decomposition

svd.Y <- svd(Y)
lambda <- svd.Y$d

# Test statistic

T <- {lambda^2/rev(cumsum(rev(lambda^2)))} [1:(M - 1)]


# Sampling

T.boot <- matrix(NA, nrow = N.Boot, ncol = M - 1)

for(Component in 1:(M - 1)){
  K <- Component - 1
  for(Boot in 1:N.Boot){
    E.boot <- matrix(rnorm((n - 1 - K)*(p - K)), nrow = n - 1 - K, ncol = p - K)
    lambda.boot <- svd(E.boot)$d
    T.boot[Boot, Component] <- lambda.boot[1]^2 / sum(lambda.boot^2)
  }
  if (mean(T.boot[, Component] > T[Component]) > alpha){break}
}

# Results
p_value <- colMeans(T.boot > matrix(rep(T, N.Boot), nrow = N.Boot, byrow = TRUE))
p_value <- na.omit(p_value)
m <- length(p_value)
result <- data.frame(matrix(c((lambda^2)[1:m], T[1:m], p_value), ncol=3))
colnames(result) <- c("lambda^2", "T", "p-value")
rownames(result) <- 1:nrow(result)
result