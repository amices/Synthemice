
library(Rcpp)
library(RcppArmadillo)
library(mice)
library(tidyverse)
library(magrittr)

set.seed(123)                   # seed
nsim <- 100                     # number os simulations
n  <- c(25, 100, 400)  # sample size
r2 <- c(0.05, 0.25)             # proportion explained variance to be divided
ratio1 <- c(1,1,1,1,1)          # over 5 predictors with (1) equal size, and
ratio2 <- c(1,2,3,4,5)          # (2) tapering size
rho1 <- diag(5)               # with either uncorrelated and correlated (r = 0.3) predictors
rho2 <- matrix(c(1, 0.3, 0.3, 0.3, 0.3,
                   0.3, 1, 0.3, 0.3, 0.3,
                   0.3, 0.3, 1, 0.3, 0.3,
                   0.3, 0.3, 0.3, 1, 0.3,
                   0.3, 0.3, 0.3, 0.3, 1), nrow = 5) 
model <- c("normal", "logit")

conditions <- expand_grid(nsim = 1:nsim, 
                          n = n,
                          r2 = r2, 
                          ratio = list(ratio1, ratio2), 
                          rho = list(rho1, rho2),
                          model = model)

make.coefficients <- function(r2, ratio, rho, model = c("normal", "logit")) {
  
  if (model == "normal") {               # given a multivariate normal model with              
    var_xb <- r2                         # standardized predictors, r2 = var(XB) / var(y)
  }
  else if (model == "logit") {           # given a logistic model with multivariate normal
    var_xb <- (r2 * pi^2 / 3) / (1 - r2) # predictors and binary outcome r2 = var(XB) / (var(XB) + pi^2/3)
  }
  sqrt(var_xb / sum(ratio %*% t(ratio) * rho)) * ratio # population coefficients of the predictors
}

x <- conditions %>%
  mutate(betas = pmap(., function(r2, ratio, rho, model, nsim, n) {
                      make.coefficients(r2, ratio, rho, model)})) %>%
  mutate(data  = pmap(., function(r2, betas, rho, n, model, nsim, ratio) {
                      generate.data(r2, betas, rho, n, model)}))
x$coefs[[301]]

make.coefficients(0.05, ratio1, rho1, "normal")
sourceCpp("~/Documents/Master_Thesis/MSBBSS/simulations/mvn_c.cpp")

object.size(x)/1000000
