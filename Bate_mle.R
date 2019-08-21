library("EnvStats")

betamle <- function(mysd){
  data <- 0.95 + rnorm(n, 0, mysd)
  data[which(data > 1)] = 1
  boots_alpha <- rep(0, m)
  boots_beta <- rep(0, m)
  for (i in 1 : m) {
    boots_data <- sample(data, replace = T)
    boots_result <- ebeta(boots_data, method = "mle") #or method = "mme"
    boots_alpha[i] <- boots_result$parameter[1]
    boots_beta[i] <- boots_result$parameter[2]
  }
  quan_alpha <- quantile(boots_alpha, c(0.025, 0.975))
  mean_alpha <- mean(boots_alpha)
  quan_beta <- quantile(boots_beta, c(0.025, 0.975))
  mean_beta <- mean(boots_beta)
  
  return(list(quan_alpha = quan_alpha, quan_beta = quan_beta, mean_alpha = mean_alpha, mean_beta = mean_beta))
}

n <- 100
m <- 10000
mysd <- seq(from = 0.001, to = 0.1, by = 0.001)
result_alpha_lo <- rep(0, length(mysd))
result_alpha_up <- rep(0, length(mysd))
result_beta_lo <- rep(0, length(mysd))
result_beta_up <- rep(0, length(mysd))
for (i in 1 : length(mysd)) {
  tmp_result <- betamle(mysd[i])
  result_alpha_lo[i] <- tmp_result$quan_alpha[1]
  result_alpha_up[i] <- tmp_result$quan_alpha[2]
  result_beta_lo[i] <- tmp_result$quan_beta[1]
  result_beta_up[i] <- tmp_result$quan_beta[2]
}
plot(mysd, result_alpha_lo, type = "l", col = "blue",xlab="(a)", ylab = "alpha_hat", ylim = c(0, 50), main = " ")
lines(mysd, result_alpha_up, col = "blue")
abline(a=15, b=0,lty=4)
abline(a=12, b=0,lty=3)
legend("topright",c("alpha_hat=15","alpha_hat=12","95% CI for alpha_hat"),cex=0.8,lty=c(4,3,1),col=c(1,1,4))

plot(mysd, result_beta_lo, type = "l", col = "blue",xlab="(b)", ylab = "beta_hat", ylim = c(0, 50), main = " ")
lines(mysd, result_beta_up, col = "blue")
abline(a=1, b=0,lty=3)
abline(a=1, b=1.5,lty=4)
legend("topright",c("beta_hat=1.5","beta_hat=1","95% CI for beta_hat"),cex=0.8,lty=c(4,3,1),col=c(1,1,4))

