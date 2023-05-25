#
require(agridat)
require(reshape2)
data(kang.peanut)
Y0 <- acast(kang.peanut, gen~env, value.var="yield", fun.aggregate = mean)
# Specify the significance level and the number of bootstrap samples

alpha <- 0.05

# Program start (no changes are needed below)

n <- nrow(Y0) ; n   
p <- ncol(Y0) ; p
M <- min(n - 1, p)

# Variables are centred, but not scaled

Y <- scale(Y0, center = TRUE, scale=TRUE)


# Singular value decomposition

svd.Y <- svd(Y)
lambda <- svd.Y$d

mu_np <- (sqrt(n-1) + sqrt(p))^2
sigma_np <- (sqrt(n-1) + sqrt(p))*(1/sqrt(n-1) + 1/sqrt(p))^(1/3)
tau_1 <- max(lambda)

l_1 <- max(svd((Y%*%diag(sqrt(rchisq(p, n))))/(n-1))$d)
W <- (l_1^2-mu_np)/sigma_np
V <- (tau_1^2-mu_np)/sigma_np
# p-value by Johnstone 2001
1-ptw(V, 1)
1-ptw(W, 1)