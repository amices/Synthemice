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

#########################################################################################
## Function for synthesizing the data, that handles partitioned data and the           ##
## complete data                                                                       ##
#########################################################################################

synthesize_2 <- function(data, parts = NULL,    # Use a function to synthesize and analyze
                       partition = FALSE, n.parts = 5,
                       method = NULL, m = 1) {   # the (synthesized) data.
  
  if (partition) {
    # If the data must be partitioned every iterations, the function creates the partitions
    partitions <- rep(1/n.parts, n.parts)
    names(partitions) <- paste0("d", 1:n.parts)
    parts <- resample_partition(data, partitions)
  }
  
  if (!is.null(parts)) {
    # If there are partitions included, synthesize synthesize the data of all partitions,
    # and rowbind the synthesized partitions
    syn_dat <- map_dfr(parts, ~ syn(.[[1]], print.flag = F, m = m, method = method)$syn) 
  }                                           
  else {
    # If it concerns the complete data as a whole, simply synthesize the complete data, 
    # and extract this data
    syn_dat <- syn(data, print.flag = F, m = m, method = method)
  }
  return(syn_dat)
}

#########################################################################################
## Execute function with default synthesizing settings                                 ##
## complete data                                                                       ##
#########################################################################################

nsim <- 1000        # number of iterations

plan(multisession) # specify parallel processing


syn_meth <- c("cart", "cart", "cart", "~ I(wgt / (hgt/100)^2)", 
              "cart", "cart", "cart", "cart", "cart")

syns_m5_complete1 <- future_map(1:nsim, ~ synthesize_2(boyscomp, m = 5, method = syn_meth),
                               .progress = TRUE,
                               .options = future_options(seed = as.integer(123)))


cutoff <- qnorm(.975)

syn5_out <- syns_m5_complete1 %>%
  map(., function(x) lm.synds(age ~ hgt + tv, data = x)) %>%
  map(., summary) %>%
  map(., coef) %>%
  map(., function(x) rownames_to_column(as.data.frame(x), var = "Variable")) %>%
  bind_rows(, .id = "sim") %>%
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




