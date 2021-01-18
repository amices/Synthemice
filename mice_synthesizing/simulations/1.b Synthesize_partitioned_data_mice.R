
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

# make sure that wgt & hgt are imputed independent of bmi
pred <- make.predictorMatrix(truth)
pred[c("wgt", "hgt"), "bmi"] <- 0

# create bootstrap samples as data.frames
boot_boys <- truth %>% bootstrap(nsim) %$% strap %>% map(as.data.frame)

# number of partitions
n.parts <- 5

# create the partitions and name them
parts <- rep(1/n.parts, n.parts)
names(parts) <- paste0("P", 1:n.parts)

# draw randomly n.parts partitions, and impute all these partitions using 
# mice with specified parameters, and complete all these subsets in long
# format. Bind these partitions and create lists, such that for 1 to m,
# the data belonging to the same imputation round are bind together (that is,
# we create m datasets of the same dimension as the original dataset).
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

# Same as the previous, but then with 5 iterations instead of 1. 
boot_parts_5 <- boot_boys %>%
  future_map(function(x) {
    resample_partition(x, parts) %>%
      map(function(y) {
        as.data.frame(y) %>%
          mice(m = 5, 
               maxit = 5,
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
