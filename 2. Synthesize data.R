#########################################################################################
## Create synthetic versions of the data.                                              ##
## build-in in mice.                                                                   ##
#########################################################################################

#########################################################################################
## Load packages, source file wherein the data is created and set seed                 ##
#########################################################################################

library(tidyverse)                            # Load tidyverse for tibbles
library(synthpop)                             # Load synthpop for data synthesizing
library(furrr)                                # furrr for parallel mapping
library(broom)                                # broom for tidy results

source("1.a Create data - boys complete.R")   # Obtain the previously synthesized data

set.seed(123)                                 # Set the seed for reproducibility

#########################################################################################
## Function for synthesizing the data, that handles partitioned data and the           ##
## complete data                                                                       ##
#########################################################################################

synthesize <- function(data, parts = NULL,    # Use a function to synthesize and analyze
                       partition = FALSE, n.parts = 5) {   # the (synthesized) data.
  
  if (partition) {
    # If the data must be partitioned every iterations, the function creates the partitions
    partitions <- rep(1/n.parts, n.parts)
    names(partitions) <- paste0("d", 1:n.parts)
    parts <- resample_partition(data, partitions)
  }
  
  if (!is.null(parts)) {
    # If there are partitions included, synthesize synthesize the data of all partitions,
    # and rowbind the synthesized partitions
    syn_dat <- map_dfr(parts, ~ syn(.[[1]],print.flag = F)$syn) 
  }                                           
  else {
    # If it concerns the complete data as a whole, simply synthesize the complete data, 
    # and extract this data
    syn_dat <- syn(data, print.flag = F)$syn
  }
  return(syn_dat)
}

nsim <- 1000        # number of iterations
plan(multisession)  # specify parallel processing

# replicate the function 1000 times on the complete data in parallel
syns_complete <- future_map_dfr(1:nsim, ~ synthesize(boyscomp),
                                .id = "sim", .progress = TRUE)
# replicate the function 1000 times on the same partitioned data
syns_part <- future_map_dfr(1:nsim, ~ synthesize(boyscomp, parts = part),
                            .id = "sim", .progress = TRUE)
# replicate the function 1000 times on a 1000 times partitioned dataset

syns_resample_partition <- future_map_dfr(1:nsim, ~ synthesize(boyscomp, partition = TRUE),
                                          .id = "sim", .progress= TRUE)

# Specify the model to be run on all synthesized datasets
syn_model <- function(data) lm(age ~ hgt + tv, data = data)

# Combine the three result objects into one output object, group it by iteration and 
# synthesizing method (complete data, initially partitioned data, iteratively partitioned
# data), nest the data. Apply the syn model to all synthesized datasets, extract the 
# coefficients and the confidence intervals, unnest the coefficients and confidence
# intervals, and add the known true estimates, and the confidence interval coverage (that
# is, how often are the true estimates included in the estimated confidence intervals).
results <- bind_rows(syns_complete, syns_part, syns_resample_partition, .id = "method") %>%
  group_by(sim, method) %>%
  nest() %>%
  mutate(model = map(data, syn_model),
         coefs = map(model, tidy),
         confint = future_map(model, confint_tidy)) %>%
  unnest(., c(coefs, confint)) %>%
  mutate(real = tidy(syn_model(boyscomp))$estimate,
         covered = real > conf.low & real < conf.high)

# Format the results so that we can nicely display them.
syn_out <- results %>%
  dplyr::select(term, real, estimate, conf.low, conf.high, covered) %>%
  group_by(method, term) %>%
  summarise("Population value" = unique(real),
            "Qbar" = mean(estimate),
            "Lower" = mean(conf.low),
            "Upper" = mean(conf.high),
            "Coverage" = mean(covered)) %>%
  mutate(method = recode(method, `1` = "Complete data", `2` = "Single partition",
                         `3` = "Nsim partitions"))
