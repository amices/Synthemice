
library(tidyverse)
library(magrittr)
library(synthpop)
library(mice)
library(furrr)

generate <- function(n) {
  X3 <- rnorm(n)
  X2 <- rnorm(n)
  X1 <- X3 * X2 + rnorm(n)
  
  Y <- X1 + X1 * X2 + X3 + rnorm(n, 0, 2)
  bind_cols(Y = Y, X1 = X1, X2 = X2, X3 = X3)
}

nsim <- 500

plan(multisession)
options <- furrr_options(seed = as.integer(123))

out <- tibble(nsim = 1:nsim)

out <- out %>%
  mutate(data = map(nsim, ~generate(100)),
         syns = future_map(data, 
                           ~synthpop::syn(.x, 
                                          m = 5,
                                          cart.minbucket = 3,
                                          cart.cp = 1e-32,
                                          print = FALSE), 
                           .options = options),
         mice = future_map(data, 
                           ~mice(.x,
                                 maxit = 5, 
                                 method = "cart", 
                                 where = matrix(1, nrow(.x), ncol(.x)),
                                 minbucket = 3,
                                 cp = 1e-32,
                                 print = FALSE),
                           .options = options))


syn$syn %>% map(~lm(Y ~ X1, .x))

out <- out %>%
  mutate(normal_lm = map(data, ~lm(Y ~ X1*X2 + X3, .x)),
         synth_lm = map(syns, ~map(.x$syn, 
                                 function(.y) lm(Y ~ X1*X2 + X3, .y)) %>% pool3.syn),
         mice_lm = map(mice, ~complete(.x, action = "all") %>%
                              map(function(.y) lm(Y ~ X1*X2 + X3, .y)) %>%
                              pool3.syn))

normal_cov <- function(true_est, est) {
  l <- confint(est)[,1]
  u <- confint(est)[,2]
  
  tibble(true = true_est,
         syn  = coef(est),
         bias = syn - true,
         cov  = l < true & true < u)
}

syn_cov <- function(true_est, pooled) {
  l <- pooled$lower
  u <- pooled$upper
  
  tibble(true = true_est,
         syn  = pooled$est,
         bias = syn - true,
         cov  = l < true & true < u )
}

out_cov <- out %>%
  mutate(normal_cov = map(normal_lm, ~normal_cov(c(0,1,0,1,1), .x)),
         synth_cov  = map(synth_lm, ~syn_cov(c(0,1,0,1,1), .x)),
         mice_cov   = map(mice_lm, ~syn_cov(c(0,1,0,1,1), .x)))

map_dfr(out$mice, 
        ~.x %>%
          complete(action = "all") %>%
          map(~lm(Y ~ X1*X2 + X3, .x)) %>%
          pool3.syn %>%
          mutate(cov = lower < c(0,1,0,1,1) & c(0,1,0,1,1) < upper)
) %>%
  group_by(term) %>%
  summarize(across(c(est, cov), mean))

map_dfr(out$syns,
        ~.x$syn %>%
          map(~lm(Y ~ X1*X2 + X3, .x)) %>%
          pool3.syn %>%
          mutate(cov = lower < c(0,1,0,1,1) & c(0,1,0,1,1) < upper)
) %>%
  group_by(term) %>%
  summarize(across(c(est, cov), mean))
