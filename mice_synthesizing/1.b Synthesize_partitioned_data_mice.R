
library(mice) # imputations
library(tidyverse) # tidy data
library(modelr) # bootstrap
library(magrittr) # pipe
library(furrr) # parallel mapping

set.seed(123) # seed for reproducibility

#create one complete dataset that allow us to work without missings
truth <- boys %>% mice(seed = 123, m = 1, print = FALSE) %>% complete()

# shorten the complete lm model, since we want to use the same model anyhow
model <- function(data) lm(wgt ~ age + hgt, data)

# run this model
truemodel <- truth %>% model

# extract the coefficients
coefs <- coef(truemodel)

# parallel processing to increase speed
plan(multisession)

# number of iterations
nsim <- 500

# cart method, all variables are imputed by means of cart
cart <- rep("cart", ncol(truth))
names(cart) <- colnames(truth)
cart['bmi'] <- "~I(wgt / (hgt/100)^2)"

pred <- make.predictorMatrix(truth)
pred[c("wgt", "hgt"), "bmi"] <- 0

boot_boys <- truth %>% bootstrap(nsim) %$% strap %>% map(as.data.frame)

n.parts <- 5

parts <- rep(1/n.parts, n.parts)
names(parts) <- paste0("P", 1:n.parts)

boot_parts <- boot_boys %>%
  future_map(function(x) {
    resample_partition(x, parts) %>%
      map(function(y) {
        as.data.frame(y) %>%
          mice(m = 5, 
               maxit = 1,
               method = cart,
               minbucket = 3,
               cp = 1e-08,
               predictorMatrix = pred,
               where = matrix(1, nrow(.), ncol(.)),
               print = F) %>% 
          mice::complete(., action = "long")
          }
        ) %>% bind_rows %>%
      plyr::dlply(~.imp)
  }, .options = future_options(seed = as.integer(123)), .progress = T)

