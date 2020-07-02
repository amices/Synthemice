#########################################################################################
## Create synthetic versions of the data.                                              ##
## build-in in mice.                                                                   ##
#########################################################################################

#########################################################################################
## Load packages, source file wherein the data is created                              ##
#########################################################################################

library(tidyverse)                            # Load tidyverse for tibbles
library(synthpop)                             # Load synthpop for data synthesizing
library(furrr)                                # furrr for parallel mapping
library(broom)                                # broom for tidy results

source("1.a Create data - boys complete.R")   # Obtain the previously synthesized data
source("1.c Functions.R")

#########################################################################################
## Execute function with default synthesizing settings                                 ##
## complete data                                                                       ##
#########################################################################################

nsim <- 1000        # number of iterations

plan(multisession) # specify parallel processing


syn_meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", 
              "cart", "cart", "cart", "cart", "cart")

syns_m5_complete1 <- future_map(1:nsim, ~ synthesize(boyscomp, m = 5, method = syn_meth),
                               .progress = TRUE,
                               .options = future_options(seed = as.integer(123)))


cutoff <- qnorm(.975)

syn5_out <- syns_m5_complete1 %>%
  map(., function(x) lm.synds(age ~ hgt + tv, data = x)) %>%
  map(., summary) %>%
  map(., coef) %>%
  map(., function(x) rownames_to_column(as.data.frame(x), var = "Variable")) %>%
  bind_rows(., .id = "sim") %>%
  rename(., "est" = "xpct(Beta)", "se" = "xpct(se.Beta)", "z" = "xpct(z)", "p" = "Pr(>|xpct(z)|)") %>%
  mutate(lower = est - cutoff*se,
         upper = est + cutoff*se) %>%
  group_by(sim) %>%
  mutate(real = coef(lm(age ~ hgt + tv, boyscomp)),
         covered = real > lower & real < upper) %>%
  group_by(Variable) %>%
  summarise("Population value" = unique(real),
            "Qbar" = mean(est),
            "Bias" = mean(est) - unique(real),
            "Lower" = mean(lower),
            "Upper" = mean(upper),
            "Coverage" = mean(covered))




