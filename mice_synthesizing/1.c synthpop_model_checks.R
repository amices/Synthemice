library(tidyverse)
library(synthpop)
library(furrr)
library(modelr)
library(magrittr)
library(broom)

source("/Users/thomvolker/Documents/Federated_imputation/synthpop_synthesizing/1.a Create data - boys complete.R")
source("/Users/thomvolker/Documents/Federated_imputation/synthpop_synthesizing/1.c Functions.R")

set.seed(1996)

format_out <- function(coefs, boot_out = NULL, syn_out = NULL) {
  
  if (!is.null(boot_out)) {
    boot_out %>%
      mutate(true = rep(coefs, nrow(boot_out)/length(coefs)),
             covered = conf.low < true & true < conf.high,
             bias = estimate - true) %>%
      ungroup() %>% group_by(term) %>%
      summarise("true" = unique(true),
                "qbar" = mean(estimate),
                "bias" = mean(bias),
                "SEbar" = mean(std.error),
                "lower" = mean(conf.low),
                "upper" = mean(conf.high),
                "ciw" = mean(upper - lower),
                "coverage" = mean(covered)) %>% print
  }
  
  if (!is.null(syn_out)) {
    syn_out %>%
      mutate(true = rep(coefs, nrow(syn_out)/length(coefs)),
             estimate = Beta.syn,
             std.error = se.Beta.syn,
             conf.low = estimate - qnorm(.975) * std.error,
             conf.high = estimate + qnorm(.975) * std.error,
             covered = conf.low < true & true < conf.high,
             bias = estimate - true) %>%
      ungroup() %>% group_by(term) %>%
      summarise("true" = unique(true),
                "qbar" = mean(estimate),
                "bias" = mean(bias),
                "SEbar" = mean(std.error),
                "lower" = mean(conf.low),
                "upper" = mean(conf.high),
                "ciw" = mean(upper - lower),
                "coverage" = mean(covered)) %>% print
  }
}

plan(multisession)

seed1 <- as.integer(123)
seed2 <- as.integer(12345)
seed3 <- as.integer(54321)

nsim <- 2000

bootstrap_samples <- modelr::bootstrap(boyscomp, nsim)

real_coefs <- coef(lm(wgt ~ age + hgt, boyscomp))

true_models <- bootstrap_samples$strap %>%
  map(function(x) as.data.frame(x) %$% lm(wgt ~ age + hgt))

boot_out <- true_models %>%
  map_dfr(function(x) tidy(x, conf.int = TRUE))

meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", "cart", "cart", "cart", "cart", "cart")

syn_models <- bootstrap_samples$strap %>%
  map(function(x) as.data.frame(x)) %>%
  future_map(function(x) {
      synthpop::syn(x, method = meth, m = 5, print.flag = FALSE) %>%
      lm.synds(wgt ~ age + hgt, .)
  }, .progress = TRUE, .options = future_options(seed = seed1))

syn_models2 <- bootstrap_samples$strap %>%
  map(function(x) as.data.frame(x)) %>%
  future_map(function(x) {
    synthpop::syn(x, method = meth, m = 5, print.flag = FALSE) %>%
      lm.synds(wgt ~ age + hgt, .)
  }, .progress = TRUE, .options = future_options(seed = seed2))

syn_models3 <- bootstrap_samples$strap %>%
  map(function(x) as.data.frame(x)) %>%
  future_map(function(x) {
    synthpop::syn(x, method = meth, m = 5, print.flag = FALSE) %>%
      lm.synds(wgt ~ age + hgt, .)
  }, .progress = TRUE, .options = future_options(seed = seed3))

syn_models_vs <- bootstrap_samples$strap %>%
  map(function(x) as.data.frame(x)) %>%
  future_map(function(x) {
    synthpop::syn(x, method = meth, m = 5, print.flag = FALSE,
                  visit.sequence = c("hgt", "wgt", "age", "bmi", 
                                     "hc", "gen", "phb", "tv", "reg")) %>%
      lm.synds(wgt ~ age + hgt, .)
  }, .progress = TRUE, .options = future_options(seed = seed1))

syn_out <- syn_models %>%
  map(function(x) summary(x, population.inference = TRUE)) %>%
  map_dfr(function(x) {
    coef(x) %>% 
      as.data.frame %>%
      rownames_to_column(var = "term")})

syn_out3 <- syn_models3 %>%
  map(function(x) summary(x, population.inference = TRUE)) %>%
  map_dfr(function(x) {
    coef(x) %>% 
      as.data.frame %>%
      rownames_to_column(var = "term")})
  

sum_out %>% 
  coef %>% 
  as.data.frame %>% 
  rownames_to_column(var = "term") %>%
  mutate(est)

coef(sum_out)
