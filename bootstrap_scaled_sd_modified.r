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


# Variables are centred, and divided by standard deviations

Y <- scale(Y0, center = TRUE, scale = TRUE)


# Singular value decomposition

svd.Y <- svd(Y)
lambda <- svd.Y$d


# Test statistic

T <- {lambda^2/rev(cumsum(rev(lambda^2)))} [1:(M - 1)]


# Sampling

T.boot <- matrix(NA, nrow = N.Boot, ncol = M - 1)
Theta.hat <- matrix(NA, ncol = 1, nrow = M - 1)
Matrix.0 <- matrix(0, nrow = n, ncol = p)
Matrix.1 <- matrix(1, nrow = n, ncol = p)

for(Component in 1:(M - 1)){
  K <- Component - 1
  Matrix.Component <- matrix(Component, nrow = n, ncol = p)
  Theta.hat <- ifelse(Matrix.Component == Matrix.1, Matrix.0, svd.Y$u[, 1:K]%*%diag(svd.Y$d[1:K], nrow = K, ncol = K)%*%t(svd.Y$v[, 1:K]))
  Variance.hat <- t(svd.Y$d[Component:M])%*%svd.Y$d[Component:M]/((n - 1 - K)*(p - K))
  for(Boot in 1:N.Boot){
    Y0.boot <- matrix(rnorm(n*p, as.vector(Theta.hat), sd = sqrt(Variance.hat)), nrow = n, ncol = p)
    Y.boot <- scale(Y0.boot, center = TRUE, scale = TRUE)
    lambda.boot <- svd(Y.boot)$d
    T.boot[Boot, Component] <- lambda.boot[Component]^2 /t(lambda.boot[Component:M])%*%lambda.boot[Component:M]
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