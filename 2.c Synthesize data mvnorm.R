library(furrr)
library(synthpop)
library(tidyverse)

bootstrap_syn <- function(data, method = NULL, formula, coefs) {
  d   <- data                                             # specify the data
  s1  <- syn(d, m = 1, method = method, print.flag = F)   # 1 synthesis
  s5  <- syn(d, m = 5, method = method, print.flag = F)   # 5 syntheses
  s10 <- syn(d, m = 10, method = method, print.flag = F)  # 10 syntheses
  
  f <- as.formula(formula)                                # the formula
  
  # Fit the lm models on all four datasets (the real data and the synthetic data)
  r <- list(Real_Sample = lm(f,d), Syn1 = lm.synds(f, s1), Syn5 = lm.synds(f, s5), Syn10 = lm.synds(f, s10))
  
  # Collect the output
  out <- map(r, function(x) as.data.frame(coef(summary(x)))) %>%
    map(., function(x) {colnames(x) <- c("Est", "SE", "stat", "P"); return(as.data.frame(x))}) %>%
    map(., function(x) rownames_to_column(x, var = "Variable")) %>%
    bind_rows(., .id = "Method") %>%
    group_by(Method) %>% 
    ungroup() %>%
    dplyr::select(., Method, Variable, Est, SE, stat, P) %>%
    mutate(Lower = Est - qnorm(.975)*SE,
           Upper = Est + qnorm(.975)*SE)
  # Return the output
  return(as.data.frame(out))
}

# Generate multicariate normal data
normal <- function(r2, ratio_beta, rho, n = 10000) {
  X <- MASS::mvrnorm(n = n, mu = rep(0, length(ratio_beta)), Sigma = rho)
  colnames(X) <- paste0("IV", 1:length(ratio_beta))
  coefs <- matrix(NA, nrow = length(ratio_beta), ncol = length(ratio_beta))
  for (i in 1:length(ratio_beta)) {
    for (j in 1:length(ratio_beta)) {
      coefs[i,j] <- ifelse(i > j, ratio_beta[i] * ratio_beta[j], 0)
    }
  }
  b <- sqrt((r2 / (sum(ratio_beta^2) + 2 * sum(coefs * rho))))
  y <- X %*% (b*ratio_beta) + rnorm(n, 0, sqrt(1 - r2))
  return(list(dat = data.frame(DV = y, X), coefs = b*ratio_beta))
}

# The ratio of the effect of the regression coefficients in the population
ratio <- c(0,1,2)
# Specify the correlation between the predictors
rho0 <- diag(length(ratio))
# Explicitly say that all off-diagonal elements are equal to zero (uncorrelated predictors)
rho0[rho0!=1] <- 0
# Test whether the function works, with an R2 of .5, the ratio of the regression coefficients
# and the correlation as specified above, and a sample size of n = 100
test_rho0 <- normal(.5, ratio_beta = ratio, rho = rho0, n = 100)
# Test whether the population regression coefficients are in line with the specified r2 (they are)
sum(test_rho0$coefs^2) + 2*sum(test_rho0$coefs[2]*test_rho0$coefs[3]*.9)

real_coefs_rho0 <- test_rho0$coefs

plan(multisession)

# method for normal synthesis (based on a linear regression model)
synds_norm <- c("norm", "norm", "norm", "norm")
# method for synthesis based on classification and regression trees
synds_cart <- c("cart", "cart", "cart", "cart")

# Generate random normal data, and run the synthesize function on it (bootstrap is a somewhat 
# inappropriate name, because we sample directly from the population). 
out_norm_rho0 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat; bootstrap_syn(data, method = synds_norm, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                .id = "sim", .progress = TRUE, .options = future_options(seed = as.integer(123)))
out_cart_rho0 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat; bootstrap_syn(data, method = synds_cart, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                .id = "sim", .progress = TRUE, .options = future_options(seed = as.integer(123)))

# Collect the results
out_norm_rho0 %>%
  mutate(RealEst = rep(real_coefs_rho0, 2000),
         Covered = Lower < real_coefs_rho0 & Upper > real_coefs_rho0,
         Bias = Est - RealEst) %>%
  ungroup() %>% group_by(Method, Variable) %>%
  summarise("Population estimate" = unique(RealEst),
            "Qbar" = mean(Est),
            "MaxSE" = max(SE),
            "MinSE" = min(SE),
            "MeanSE" = mean(SE),
            "Bias" = mean(Bias),
            "Lower" = mean(Lower),
            "Upper" = mean(Upper),
            "Coverage" = mean(Covered)) %>%
  arrange(factor(Method, levels = c("Real_Sample", "Syn1", "Syn5", "Syn10"))) -> summary_norm_rho0

out_cart_rho0 %>%
  mutate(RealEst = rep(real_coefs_rho0, 2000),
         Covered = Lower < real_coefs_rho0 & Upper > real_coefs_rho0,
         Bias = Est - RealEst) %>%
  ungroup() %>% group_by(Method, Variable) %>%
  summarise("Population estimate" = unique(RealEst),
            "Qbar" = mean(Est),
            "MaxSE" = max(SE),
            "MinSE" = min(SE),
            "MeanSE" = mean(SE),
            "Bias" = mean(Bias),
            "Lower" = mean(Lower),
            "Upper" = mean(Upper),
            "Coverage" = mean(Covered)) %>%
  arrange(factor(Method, levels = c("Real_Sample", "Syn1", "Syn5", "Syn10"))) -> summary_cart_rho0

ratio <- c(0,1,2)
rho50 <- diag(length(ratio))
rho50[rho50!=1] <- .5        # Now we use a correlation between the predictors of .5
test_rho50 <- normal(.5, ratio_beta = ratio, rho = rho50, n = 100)

sum(test_rho50$coefs^2) + 2*sum(test_rho50$coefs[2]*test_rho50$coefs[3]*.9)

real_coefs_rho50 <- test_rho50$coefs

# And run the same function as above, but now with correlated predictors
out_norm_rho50 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat; bootstrap_syn(data, method = synds_norm, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                 .id = "sim", .progress = TRUE, .options = future_options(seed = as.integer(123)))
out_cart_rho50 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat; bootstrap_syn(data, method = synds_cart, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                 .id = "sim", .progress = TRUE, .options = future_options(seed = as.integer(123)))
# And collect the output
out_norm_rho50 %>%
  mutate(RealEst = rep(real_coefs_rho50, 2000),
         Covered = Lower < real_coefs_rho50 & Upper > real_coefs_rho50,
         Bias = Est - RealEst) %>%
  ungroup() %>% group_by(Method, Variable) %>%
  summarise("Population estimate" = unique(RealEst),
            "Qbar" = mean(Est),
            "MaxSE" = max(SE),
            "MinSE" = min(SE),
            "MeanSE" = mean(SE),
            "Bias" = mean(Bias),
            "Lower" = mean(Lower),
            "Upper" = mean(Upper),
            "Coverage" = mean(Covered)) %>%
  arrange(factor(Method, levels = c("Real_Sample", "Syn1", "Syn5", "Syn10"))) -> summary_norm_rho50

out_cart_rho50 %>%
  mutate(RealEst = rep(real_coefs_rho50, 2000),
         Covered = Lower < real_coefs_rho50 & Upper > real_coefs_rho50,
         Bias = Est - RealEst) %>%
  ungroup() %>% group_by(Method, Variable) %>%
  summarise("Population estimate" = unique(RealEst),
            "Qbar" = mean(Est),
            "MaxSE" = max(SE),
            "MinSE" = min(SE),
            "MeanSE" = mean(SE),
            "Bias" = mean(Bias),
            "Lower" = mean(Lower),
            "Upper" = mean(Upper),
            "Coverage" = mean(Covered)) %>%
  arrange(factor(Method, levels = c("Real_Sample", "Syn1", "Syn5", "Syn10"))) -> summary_cart_rho50

## Combine the output
synth_by_lm <- bind_rows("Rho = 0" = summary_norm_rho0, "Rho = .5" = summary_norm_rho50, .id = "Correlation predictors")
synth_by_cart <- bind_rows("Rho = 0" = summary_cart_rho0, "Rho = .5" = summary_cart_rho50, .id = "Correlation predictors")





