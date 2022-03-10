library(ggplot2)
library(openxlsx)
library(mvtnorm)

dane <- read.xlsx("Dane_wyplaty.xlsx", sheet=1, startRow = 1, colNames = TRUE, rowNames = FALSE,
                  detectDates = FALSE)


dane <- read.csv("weather.csv")
y <- as.matrix(cbind("median_house_value" = dane[,7]))
X <- as.matrix(cbind("const"=matrix(1,nrow = nrow(y), 1), dane[,c(-4,-6,-7,-8)]))
tytul_y <- colnames(dane[7])
tytuly <- colnames(X)
tytuly <- append(tytuly, "tau")

BMNRL <- function(y, X, burnin, mcmc, a, C, n0, s0, plots){
  T <- nrow(y)
  k <- ncol(X)
  C <- diag(k)*C
  beta0 <- solve(t(X)%*%X)%*%t(X)%*%y
  tau0 <- 1
  N_total <- burnin+mcmc
  beta <- beta0
  tau <- tau0
  S_beta <- t(y-X%*%beta0)%*%(y-X%*%beta0)
  beta_matrix <- matrix(0, nrow = N_total, ncol = k)
  tau_matrix <- matrix(0, nrow = N_total, ncol = 1)
  for(i in 1:N_total){
    n_fala <- n0+T
    s_fala <- s0+S_beta[i]
    tau <- rgamma(1, shape=n_fala/2, rate=s_fala/2)
    C_fala <- C+tau*t(X)%*%X
    a_fala <- solve(C_fala)%*%(C%*%a+tau*t(X)%*%X%*%beta0)
    beta <- t(rmvnorm(1, mean=a_fala, sigma = solve(C_fala)))
    S_beta[i+1] <- t(y-X%*%beta)%*%(y-X%*%beta)
    beta_matrix[i,] <- beta
    tau_matrix[i] <- tau
  }
  beta_matrix <- cbind(beta_matrix, tau_matrix)
  colnames(beta_matrix) <- tytuly
  #histogramy a post + gêstoœæ a priori
  if(plots == 1){
    for(i in 1:(k+1)){
      hist1 <- hist(beta_matrix[,i], probability=TRUE, main = tytuly[i])
      lines(hist1$mids, dnorm(hist1$mids, mean = 0, sd = (n0)^(-0.5)))
      dev.print(pdf,paste("Histogram",tytuly[i],".pdf"))
      dev.off()
    }
    #analiza zbieznosci
    #trace plots
    for(i in 1:(k+1)){
      plot.ts(beta_matrix[,i], main = tytuly[i])
      dev.print(pdf,paste("Trace plot",tytuly[i],".pdf"))
      dev.off()
    }
    #ACF
    for(i in 1:(k+1)){
      acf(beta_matrix[,i], main = tytuly[i])
      dev.print(pdf,paste("ACF",tytuly[i],".pdf"))
      dev.off()
    }
    #srednie ergodyczne
    for(i in 1:(k+1)){
      c1 <- cumsum(beta_matrix[,i])/seq_along(beta_matrix[,i])
      plot(c1, type="l", panel.first=abline(h=mean(beta_matrix[,i]), col="red"), las=1, main = tytuly[i])
      dev.print(pdf,paste("Œrednie ergodyczne",tytuly[i],".pdf"))
      dev.off()
    }
    #cumsumy
    #colnames(beta_matrix) <- tytuly
    cumuplot(beta_matrix, probs=c(0.5))
    dev.print(pdf,paste("cumuplot",tytuly[i],".pdf"))
    dev.off()
  }
  
  #charakterystyki a post
  expected_values <- apply(beta_matrix, 2, mean)
  medians <- apply(beta_matrix, 2, median)
  sd <- apply(beta_matrix, 2, sd)
  modalnas <- function(x) {
    temp <- unique(x)
    temp[which.max(tabulate(match(x, temp)))]
  }
  mod <- apply(beta_matrix, 2, modalnas)
  
  ci_95 <- apply(beta_matrix, 2, quantile, probs = c(0.025, 0.975))
  
  HPD <- HPDinterval(as.mcmc(beta_matrix), prob=0.95)
  #porownywanie modeli
  LogLike <- matrix(0, nrow = ncol(X)-1, 1)
  for(i in 2:ncol(X)){
    post <- MCMCregress(y~X[,c(-i)], data=dane, burnin=burnin,
                        mcmc=mcmc, b0=a, B0=C, c0=n0, d0=s0,
                        marginal.likelihood="Chib95")
    LogLike[i] <- attr(post, "logmarglike")
  }
  
  
  return_list <- list(Betas = beta_matrix, S_beta = S_beta, tau = tau_matrix,
                      Ex = expected_values, Mediana = medians, SD = sd, Mod = mod, ci_95 = ci_95, HPD = HPD,
                      LogLike = LogLike)
}

posterior <- BMNRL(y = y, X = X, burnin = 1000, mcmc=10000, a = matrix(0, nrow = ncol(X), ncol = 1), 
                   C = 0.01, n0 = 0.001, s0 = 0.001, plots = 0)


posterior2 <- MCMCregress(Wyplaty~trend+Q1+Q2+Q3+Bankomaty+Karty, data=dane, burnin=10000,
                          mcmc=50000, verbose = 500, b0=0, B0=0.01, c0=0.001, d0=0.001,
                          marginal.likelihood="Chib95")

for(i in 2:ncol(X)){
  post <- MCMCregress(y~X, data=dane, burnin=burnin,
                      mcmc=mcmc, b0=a, B0=C, c0=n0, d0=s0,
                      marginal.likelihood="Chib95")
  LogLike[i] <- attr(post, "logmarglike")
}


post <- MCMCregress(tytul_y~., data=dane, burnin=burnin,
                    mcmc=mcmc, b0=a, B0=C, c0=n0, d0=s0,
                    marginal.likelihood="Chib95")

