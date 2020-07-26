library(mice)
library(tidyverse)
library(modelr)
library(magrittr)
library(furrr)

set.seed(123)

truth <- boys %>% mice(seed = 123, m = 1, print = FALSE) %>% complete()

model <- function(data) lm(wgt ~ age + hgt, data)

truemodel <- truth %>% model

coefs <- coef(truemodel)

## Approach one - overwrite all of the data

plan(multisession)

nsim <- 10

def <- rep("pmm", ncol(truth))
names(def) <- colnames(truth)
def['bmi'] <- "~I(wgt / (hgt/100)^2)"
def[c('gen', 'phb')] <- "polr"
def['reg'] <- "polyreg"


cart <- rep("cart", ncol(truth))
names(cart) <- colnames(truth)

pred <- make.predictorMatrix(truth)
pred[c("wgt", "hgt"), "bmi"] <- 0

syns_def <- future_map(1:nsim, ~ {
  truth %>% mice(m = 5, 
                 method = def, 
                 predictorMatrix = pred, 
                 where = matrix(TRUE, nrow(truth), ncol(truth)), 
                 print = F)
}, .options = future_options(seed = as.integer(123)), .progress = TRUE, .id = "syn")



syns_cart <- future_map(1:nsim, ~ {
  truth %>% mice(m = 5, 
                 method = cart,
                 predictorMatrix = pred,
                 where = matrix(TRUE, nrow(truth), ncol(truth)),
                 print = F)
}, .options = future_options(seed = as.integer(123)), .progress = TRUE, .id = "syn")

cart_wgt_hgt <- future_map(1:nsim, ~ {
  truth %>% mice(m = 5,
                 method = cart,
                 predictorMatrix = pred,
                 where = matrix(TRUE, nrow(truth), ncol(truth)),
                 visitSequence = c('age', 'wgt', 'hgt', 'bmi', 'hc', 'gen', 'phb', 'tv', 'reg'),
                 print = F)
}, .options = future_options(seed = as.integer(123)), .progress = T, .id = "syn")


bootstrap_boys <- bootstrap(truth, nsim) %$% strap %>% map(as.data.frame)

boot_cart <- bootstrap_boys %>% 
  future_map(function(x) {
    mice(x, 
         m = 5,
         method = cart,
         predictorMatrix = pred,
         where = matrix(TRUE, nrow(truth), ncol(truth)),
         print = F)
    }, .options = future_options(seed = as.integer(123)), .progress = T)
