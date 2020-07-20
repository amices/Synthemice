

library(tidyverse)
library(synthpop)
library(furrr)

source("1.a Create data - boys complete.R")
source("1.c Functions.R")

plan(multisession)

seed <- as.integer(12345)

# Specify a small function for inside the future_map commands
synth_bootstrapped_data <- function(p, data, formula, pop.inf = T, method = NULL, M = 1, visit = (1:ncol(data))) {
  part <- data %>% resample_bootstrap %>% data.frame %>% resample_partition(p)
  syn_part_data(partitioned_data = part, M = M, formula = formula, pop.inf = pop.inf, method = method, visit = visit)
}

# Specify the number of partitions (5 equally large partitions)
n.parts <- 5
p <- rep(1 / n.parts, n.parts)
names(p) <- paste0("p", 1:n.parts)

nsim <- 1000

f <- wgt ~ age + hgt

meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", "cart", "cart", "cart", "cart", "cart")

wgt_hgt_swapped <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, pop.inf = T, M = 10,
                                                                         method = meth,
                                                                         visit = c("age", "wgt", "hgt", "bmi", "hc", "gen", "phb", "tv", "reg")),
                                       .id = "sim", .progress = T, .options = future_options(seed = seed))


get_coefs <- function(data) {
  data %>% 
    filter(Method == "Real_Sample") %>% 
    group_by(Variable) %>%
    summarise(coef = mean(Est)) %>%
    pull(coef, name = Variable)}


summary_wgt_hgt <- wgt_hgt_swapped %>% print_results(., get_coefs(.))




