source("1.a Create data - boys complete.R")
source("1.c Functions.R")

library(furrr)
library(tidyverse)
library(synthpop)

plan(multisession)

nsim <- 1000

r2 <- .5
ratio <- c(0,1,2)
rho <- diag(length(ratio))
n1 <- 100
real_coefs <- normal(r2,ratio,rho,n1)$coefs

out_norm <- future_map_dfr(1:nsim, function(x) {data <- normal(r2, ratio, rho, n1)$dat
                                                normal_syn(data, formula = DV ~ -1 + IV1 + IV2 + IV3, pop.inf = T)},
                           .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

summary_out_norm <- out_norm %>% filter(Method == "Real_Sample" | Method == "Syn5") %>%
  print_results(., real_coefs)

p <- c(d1 = .2, d2 = .2, d3 = .2, d4 = .2, d5 = .2)

out_norm_part_n100 <- future_map_dfr(1:nsim, function(x) {data <- resample_partition(normal(r2, ratio, rho, n1)$dat, p)
                                                          syn_part_data(data, M = 5, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

summary_out_norm_part_n100 <- print_results(out_norm_part_n100, coefs = real_coefs)

n2 <- 1000

out_norm_part_n1000 <- future_map_dfr(1:nsim, function(x) {data <- resample_partition(normal(r2, ratio, rho, n2)$dat, p)
                                                           syn_part_data(data, M = 5, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                      .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

summary_out_norm_part_n1000 <- print_results(out_norm_part_n1000, coefs = real_coefs)

n3 <- 10000

out_norm_part_n10000 <- future_map_dfr(1:nsim, function(x) {data <- resample_partition(normal(r2, ratio, rho, n3)$dat, p)
                                                            syn_part_data(data, M = 5, formula = DV ~ -1 + IV1 + IV2 + IV3)},
                                       .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

summary_out_norm_part_n10000 <- print_results(out_norm_part_n10000, coefs = real_coefs)

out_same_part <- future_map_dfr(1:nsim, function(x) syn_part_data(part, M = 5, formula = wgt ~ hgt + age),
                      .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

out_nsim_parts <- future_map_dfr(1:nsim, function(x) {part <- resample_partition(boyscomp, partitions)
                                                      syn_part_data(part, M = 5, formula = wgt ~ hgt + age)},
                                 .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

boyscoefs <- coef(lm(wgt ~ hgt + age, data = boyscomp))

summary_out_same_part <-  print_results(out_same_part, coefs = boyscoefs)
summary_out_nsim_parts <- print_results(out_nsim_parts, coefs = boyscoefs)

