

library(tidyverse)
library(synthpop)
library(furrr)

source("1.a Create data - boys complete.R")
source("1.c Functions.R")

plan(multisession)

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

boys_cart_results_10 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, pop.inf = T, M = 10),
                                    .id = "sim", .progress = T, .options = future_options(seed = as.integer(125)))

meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", "cart", "cart", "cart", "cart", "cart")

vs <- c(2,3,4,1,5,6,7,8,9)

boys_vs_results_10 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, method = meth, 
                                                                     pop.inf = T, M = 10, visit = vs), 
                                   .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

boys_cart_results_20 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, pop.inf = T, M = 20),
                                    .id = "sim", .progress = T, .options = future_options(seed = as.integer(125)))

boys_vs_results_20 <- future_map_dfr(1:nsim, ~ synth_bootstrapped_data(p = p, data = boyscomp, formula = f, method = meth, 
                                                                    pop.inf = T, M = 20, visit = vs), 
                                  .id = "sim", .progress = T, .options = future_options(seed = as.integer(123)))

get_coefs <- function(data) {
  data %>% 
    filter(Method == "Real_Sample") %>% 
    group_by(Variable) %>%
    summarise(coef = mean(Est)) %>%
    pull(coef, name = Variable)}

pop_coefs_cart_10 <- get_coefs(boys_cart_results_10)
pop_coefs_vs_10 <- get_coefs(boys_vs_results_10)
pop_coefs_cart_20 <- get_coefs(boys_cart_results_20)
pop_coefs_vs_20 <- get_coefs(boys_vs_results_20)

summary_cart_10 <- print_results(boys_cart_results_10, pop_coefs_cart_10)
summary_vs_10 <- print_results(boys_vs_results_10, pop_coefs_vs_10)
summary_cart_20 <- print_results(boys_cart_results_20, pop_coefs_cart_20)
summary_vs_20 <- print_results(boys_vs_results_20, pop_coefs_vs_20)

# 
# 
# 
# 
# boys_cart_results %>% 
#   filter(Variable == "age" & Method == "Syn") %>%
#   ggplot(mapping = aes(x = Est, color = Method, fill = Method)) +
#   geom_histogram(alpha = .4) + 
#   scale_color_brewer(palette = "Set1") +
#   scale_fill_brewer(palette = "Set1") +
#   ggtitle("Age") +
#   theme_classic()
# 
# 
# 
# 
# unique(boys_cart_results$Variable)
