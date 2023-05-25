rm(list = ls())
set.seed(42)
require(agridat)
require(reshape2)
data(kang.peanut)
Y0 <- acast(kang.peanut, gen~env, value.var="yield", fun.aggregate = mean)
yield_t_statistic <- function(Y0, method, B) {
    n <- nrow(Y0)
    p <- ncol(Y0)
    M <- min(n - 1, p)
    T.boot <- matrix(NA, nrow = B, ncol = M - 1)
    Matrix.0 <- matrix(0, nrow = n, ncol = p)
    Matrix.1 <- matrix(1, nrow = n, ncol = p)
    if (method == "scale-sd"){
        Y <- scale(Y0, center = TRUE, scale = TRUE)
        svd.Y <- svd(Y)
        lambda <- svd.Y$d
        for(Component in 1:(M - 1)) {
            K <- Component -1
            Matrix.Component <- matrix(Component, nrow = n, ncol = p)
            Theta.hat <- ifelse(Matrix.Component == Matrix.1, Matrix.0, svd.Y$u[, 1:K]%*%diag(svd.Y$d[1:K], nrow = K, ncol = K)%*%t(svd.Y$v[, 1:K]))
            Variance.hat <- t(svd.Y$d[Component:M])%*%svd.Y$d[Component:M]/((n - 1 - K)*(p - K))
            for(Boot in 1:B){
                Y0.boot <- matrix(rnorm(n*p, as.vector(Theta.hat), sd = sqrt(Variance.hat)), nrow = n, ncol = p)
                Y.boot <- scale(Y0.boot, center = TRUE, scale = TRUE)
                lambda.boot <- svd(Y.boot)$d
                T.boot[Boot, Component] <- lambda.boot[Component]^2 /t(lambda.boot[Component:M])%*%lambda.boot[Component:M]
            }
        }

    }
    
    if (method == "scale-mean") {
        means.Y0 <- apply(Y0, 2, mean)
        Y <- scale(Y0, center = TRUE, scale = means.Y0)
        svd.Y <- svd(Y)
        lambda <- svd.Y$d        
        J <- matrix(rep(1, n), ncol = 1)
        A0 <- kronecker(matrix(means.Y0, nrow = 1), J)
        D0 <- diag(means.Y0)
        for(Component in 1:(M - 1)) {
            K <- Component - 1
            Matrix.Component <- matrix(Component, nrow = n, ncol = p)
            Theta.hat <- ifelse(Matrix.Component == Matrix.1, Matrix.0, svd.Y$u[, 1:K]%*%diag(svd.Y$d[1:K], nrow = K, ncol = K)%*%t(svd.Y$v[, 1:K]))
            Variance.hat <- t(svd.Y$d[Component:M])%*%svd.Y$d[Component:M]/((n - 1 - K)*(p - K))
            for(Boot in 1:B){
                R0.boot <- matrix(rnorm(n*p, as.vector(Theta.hat), sd = sqrt(Variance.hat)), nrow = n, ncol = p)
                Y0.boot <- A0 + R0.boot%*%D0
                means.Y0.boot <- apply(Y0.boot, 2, mean)
                Y.boot <- scale(Y0.boot, center = TRUE, scale = means.Y0.boot)
                lambda.boot <- svd(Y.boot)$d
                T.boot[Boot, Component] <- lambda.boot[Component]^2 /t(lambda.boot[Component:M])%*%lambda.boot[Component:M]
            }
        }
    }
    
    
    T.boot

}

simulation_test <- function(Y0, method, alpha, B.boot, B.resample) {
    t_boot_table <- yield_t_statistic(Y0, method, B.boot)
    ci <- apply(t_boot_table, 2, quantile, probs = c(alpha/2, 1-alpha/2))
    t_resample_table <- yield_t_statistic(Y0, method, B.resample)
    error <- rowMeans(apply(t_resample_table, 1, function(row) (row < ci[1,]))) + rowMeans(apply(t_resample_table, 1, function(row) (row > ci[2,])))
    result <- data.frame(error)
    rownames(result) = 0:(length(error)-1)
    colnames(result) = method
    result
}
result.mean <- simulation_test(Y0, "scale-mean", 0.05, 10000, 1000)
result.sd <- simulation_test(Y0, "scale-sd", 0.05, 10000, 1000)
cbind(result.mean, result.sd)