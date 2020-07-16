

library(tidyverse)
library(synthpop)
library(furrr)

source("1.a.2 Create data - boys7482 complete.R")
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

boys_cart_results_10_7482 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp7482, formula = f, pop.inf = T, M = 10),
                                       .id = "sim", .progress = T, .options = future_options(seed = seed))

# meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", "cart", "cart", "cart", "cart", "cart")
# 
# vs <- c(2,3,4,1,5,6,7,8,9)
# 
# boys_vs_results_10 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, method = meth, 
#                                                                        pop.inf = T, M = 10, visit = vs), 
#                                      .id = "sim", .progress = T, .options = future_options(seed = seed))
# 
# boys_cart_results_20 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, pop.inf = T, M = 20),
#                                        .id = "sim", .progress = T, .options = future_options(seed = seed))
# 
# boys_vs_results_20 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, method = meth, 
#                                                                        pop.inf = T, M = 20, visit = vs), 
#                                      .id = "sim", .progress = T, .options = future_options(seed = seed))
# 
# boys_vs2_results_10 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, method = meth, 
#                                                                         pop.inf = T, M = 10), 
#                                       .id = "sim", .progress = T, .options = future_options(seed = seed))

get_coefs <- function(data) {
  data %>% 
    filter(Method == "Real_Sample") %>% 
    group_by(Variable) %>%
    summarise(coef = mean(Est)) %>%
    pull(coef, name = Variable)}


summary_cart_10_7482 <- boys_cart_results_10_7482 %>% print_results(., get_coefs(.))


# summary_vs_10 <- boys_vs_results_10 %>% print_results(., get_coefs(.))
# summary_cart_20 <- boys_cart_results_20 %>% print_results(., get_coefs(.))
# summary_vs_20 <- boys_vs_results_20 %>% print_results(., get_coefs(.))

# boys_vs2_results_10 %>% print_results(., get_coefs(.))
# 
# future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, 
#                                                  pop.inf = T, M = 10), 
#                .id = "sim", .progress = T, .options = future_options(seed = seed)) %>% print_results(., get_coefs(.))