library(furrr)
library(synthpop)
library(tidyverse)

source("1.c Functions.R")

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
sum(test_rho0$coefs^2) + 2*sum(test_rho0$coefs[2]*test_rho0$coefs[3]*0)

real_coefs_rho0 <- test_rho0$coefs

plan(multisession)

seed <- as.integer(12345)

# method for normal synthesis (based on a linear regression model)
synds_norm <- c("norm", "norm", "norm", "norm")
# method for synthesis based on classification and regression trees
synds_cart <- c("cart", "cart", "cart", "cart")

# Generate random normal data, and run the synthesize function on it (bootstrap is a somewhat 
# inappropriate name, because we sample directly from the population). 
out_norm_rho0 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat; normal_syn(data, method = synds_norm, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                .id = "sim", .progress = TRUE, .options = future_options(seed = seed))
out_cart_rho0 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat; normal_syn(data, method = synds_cart, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                .id = "sim", .progress = TRUE, .options = future_options(seed = seed))


summary_norm_rho0 <- print_results(out_norm_rho0, real_coefs_rho0)
summary_cart_rho0 <- print_results(out_cart_rho0, real_coefs_rho0)


ratio <- c(0,1,2)
rho50 <- diag(length(ratio))
rho50[rho50!=1] <- .5        # Now we use a correlation between the predictors of .5
test_rho50 <- normal(.5, ratio_beta = ratio, rho = rho50, n = 100)

sum(test_rho50$coefs^2) + 2*sum(test_rho50$coefs[2]*test_rho50$coefs[3]*.9)

real_coefs_rho50 <- test_rho50$coefs

# And run the same function as above, but now with correlated predictors
out_norm_rho50 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat; normal_syn(data, method = synds_norm, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                 .id = "sim", .progress = TRUE, .options = future_options(seed = seed))
out_cart_rho50 <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat; normal_syn(data, method = synds_cart, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                 .id = "sim", .progress = TRUE, .options = future_options(seed = seed))
# And collect the output
summary_norm_rho50 <- print_results(out_norm_rho50, real_coefs_rho50)
summary_cart_rho50 <- print_results(out_cart_rho50, real_coefs_rho50)

## Combine the output
synth_by_lm <- bind_rows("Rho = 0" = summary_norm_rho0, "Rho = .5" = summary_norm_rho50, .id = "Correlation predictors")
synth_by_cart <- bind_rows("Rho = 0" = summary_cart_rho0, "Rho = .5" = summary_cart_rho50, .id = "Correlation predictors")



#########################################################################################
## Specify population.inference = TRUE                                                 ##
#########################################################################################


# Generate random normal data, and run the synthesize function on it (bootstrap is a somewhat 
# inappropriate name, because we sample directly from the population). 
out_norm_rho0_pop <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat 
                                                        normal_syn(data, method = synds_norm, 
                                                                      formula = DV ~ -1 + IV1 + IV2 + IV3,
                                                                      pop.inf = T)},
                                    .id = "sim", .progress = TRUE, .options = future_options(seed = seed))
out_cart_rho0_pop <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho0, n = 100)$dat
                                                        normal_syn(data, method = synds_cart, 
                                                                      formula = DV ~ -1 + IV1 + IV2 + IV3,
                                                                      pop.inf = T)},
                                    .id = "sim", .progress = TRUE, .options = future_options(seed = seed))

# Collect the results
summary_norm_rho0_pop <- print_results(out_norm_rho0_pop, real_coefs_rho0)
summary_cart_rho0_pop <- print_results(out_cart_rho0_pop, real_coefs_rho0)

# Generate random normal data, and run the synthesize function on it (bootstrap is a somewhat 
# inappropriate name, because we sample directly from the population). 
out_norm_rho50_pop <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat 
                                                         normal_syn(data, method = synds_norm, 
                                                                       formula = DV ~ -1 + IV1 + IV2 + IV3,
                                                                       pop.inf = T)},
                                     .id = "sim", .progress = TRUE, .options = future_options(seed = seed))
out_cart_rho50_pop <- future_map_dfr(1:500, function(x) {data <- normal(.5, ratio, rho50, n = 100)$dat
                                                         normal_syn(data, method = synds_cart,
                                                                       formula = DV ~ -1 + IV1 + IV2 + IV3,
                                                                       pop.inf = T)},
                                     .id = "sim", .progress = TRUE, .options = future_options(seed = seed))

summary_norm_rho50_pop <- print_results(out_norm_rho50_pop, real_coefs_rho50)
summary_cart_rho50_pop <- print_results(out_cart_rho50_pop, real_coefs_rho50)

synth_by_lm_pop <- bind_rows("Rho = 0" = summary_norm_rho0_pop, "Rho = .5" = summary_norm_rho50_pop, .id = "Correlation predictors")
synth_by_cart_pop <- bind_rows("Rho = 0" = summary_cart_rho0_pop, "Rho = .5" = summary_cart_rho50_pop, .id = "Correlation predictors")

