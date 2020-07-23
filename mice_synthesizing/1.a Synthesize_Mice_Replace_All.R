library(mice)
library(tidyverse)
library(modelr)
library(magrittr)
library(furrr)

truth <- boys %>% mice(seed = 123, m = 1, print = FALSE) %>% complete()

model <- function(data) lm(wgt ~ age + hgt, data)

truemodel <- truth %>% model

coefs <- coef(truemodel)

## Approach one - overwrite all of the data

plan(multisession)

nsim <- 200

def <- rep("pmm", ncol(truth))
names(def) <- colnames(truth)
def['bmi'] <- "~I(wgt / (hgt/100)^2)"
def[c('gen', 'phb')] <- "polr"
def['reg'] <- "polyreg"


cart <- rep("cart", ncol(truth))
names(cart) <- colnames(truth)

pred <- make.predictorMatrix(truth)
pred[c("wgt", "hgt"), "bmi"] <- 0

all_syns_def <- future_map(1:nsim, ~ {
  truth %>% mice(m = 5, 
                 method = def, 
                 predictorMatrix = pred, 
                 where = matrix(TRUE, nrow(truth), ncol(truth)), 
                 print = F)
}, .options = future_options(seed = as.integer(123)), .progress = TRUE, .id = "syn")



all_syns_cart <- future_map(1:nsim, ~ {
  truth %>% mice(m = 5, 
                 method = cart,
                 predictorMatrix = pred,
                 where = matrix(TRUE, nrow(truth), ncol(truth)),
                 print = F)
}, .options = future_options(seed = as.integer(123)), .progress = TRUE, .id = "syn")


est_def <- all_syns_def %>% map(function(x) with(x, lm(wgt ~ age + hgt)) %>% pool %>% summary)

est_cart <- all_syns_cart %>% map(function(x) with(x, lm(wgt ~ age + hgt)) %>% pool %>% summary)

CIs_def <- map(est_def, function(x) {
  
  var      <- x[,1]
  true_est <- coefs
  est      <- x[,2]
  true_se  <- sqrt(diag(vcov(truemodel)))
  se       <- x[,3]
  lower    <- x[,2] + x[,3] * qt(.025, x[,5])
  upper    <- x[,2] + x[,3] * qt(.975, x[,5])
  cov      <- lower < coefs & coefs < upper
  
  bind_cols(var = var, true_est = true_est, est = est, true_se = true_se, se = se, lower = lower, upper = upper, cov = cov)
})

CIs_cart <- map(est_cart, function(x) {
  
  var      <- x[,1]
  true_est <- coefs
  est      <- x[,2]
  true_se  <- sqrt(diag(vcov(truemodel)))
  se       <- x[,3]
  lower    <- x[,2] + x[,3] * qt(.025, x[,5])
  upper    <- x[,2] + x[,3] * qt(.975, x[,5])
  cov      <- lower < coefs & coefs < upper
  
  bind_cols(var = var, true_est = true_est, est = est, true_se = true_se, se = se, lower = lower, upper = upper, cov = cov)
})

results_def <- CIs_def %>% bind_rows %>% group_by(var) %>%
  summarise("True Est" = unique(true_est),
            "Syn Est"  = mean(est),
            "Bias"     = mean(est - true_est),
            "True SE"  = unique(true_se),
            "Syn SE"   = mean(se),
            "CIW"      = mean(upper - lower),
            "Coverage" = mean(cov))

results_cart <- CIs_cart %>% bind_rows %>% group_by(var) %>%
  summarise("True Est" = unique(true_est),
            "Syn Est"  = mean(est),
            "Bias"     = mean(est - true_est),
            "True SE"  = unique(true_se),
            "Syn SE"   = mean(se),
            "CIW"      = mean(upper - lower),
            "Coverage" = mean(cov))

bind_rows("Mice default" = results_def,
          "Mice with cart" = results_cart, .id = "Imputation method")



## Approach two

random_missings <- function(data, prop.missing = .5) {
  
  n <- nrow(data); k <- ncol(data)
  
  indicators <- c(rep(1, ceiling(k*prop.missing)), rep(0, floor(k*prop.missing)))
  
  replicate(n, sample(indicators)) %>% t
}

impute <- function(data, meth, pred) {
  data %>% mice(m = 5, method = meth, predictorMatrix = pred, where = random_missings(.), 
                print = F) %$% lm(wgt ~ age + hgt) %>% pool %>% summary
}

est_def_50 <- future_map(1:nsim, ~ impute(truth, meth = def, pred = pred),
                             .progress = TRUE, .options = future_options(seed = as.integer(123)))

est_cart_50 <- future_map(1:nsim, ~ impute(truth, meth = cart, pred = pred),
                             .progress = TRUE, .options = future_options(seed = as.integer(123)))


CIs_def_50 <- map(est_def_50, function(x) {
  
  var      <- x[,1]
  true_est <- coefs
  est      <- x[,2]
  true_se  <- sqrt(diag(vcov(truemodel)))
  se       <- x[,3]
  lower    <- x[,2] + x[,3] * qt(.025, x[,5])
  upper    <- x[,2] + x[,3] * qt(.975, x[,5])
  cov      <- lower < coefs & coefs < upper
  
  bind_cols(var = var, true_est = true_est, est = est, true_se = true_se, se = se, lower = lower, upper = upper, cov = cov)
})

CIs_cart_50 <- map(est_cart_50, function(x) {
  
  var      <- x[,1]
  true_est <- coefs
  est      <- x[,2]
  true_se  <- sqrt(diag(vcov(truemodel)))
  se       <- x[,3]
  lower    <- x[,2] + x[,3] * qt(.025, x[,5])
  upper    <- x[,2] + x[,3] * qt(.975, x[,5])
  cov      <- lower < coefs & coefs < upper
  
  bind_cols(var = var, true_est = true_est, est = est, true_se = true_se, se = se, lower = lower, upper = upper, cov = cov)
})

results_def_50 <- CIs_def_50 %>% bind_rows %>% group_by(var) %>%
  summarise("True Est" = unique(true_est),
            "Syn Est"  = mean(est),
            "Bias"     = mean(est - true_est),
            "True SE"  = unique(true_se),
            "Syn SE"   = mean(se),
            "CIW"      = mean(upper - lower),
            "Coverage" = mean(cov))

results_cart_50 <- CIs_cart_50 %>% bind_rows %>% group_by(var) %>%
  summarise("True Est" = unique(true_est),
            "Syn Est"  = mean(est),
            "Bias"     = mean(est - true_est),
            "True SE"  = unique(true_se),
            "Syn SE"   = mean(se),
            "CIW"      = mean(upper - lower),
            "Coverage" = mean(cov))

bind_rows("Mice default" = results_def_50,
          "Mice with cart" = results_cart_50, .id = "Imputation method")

# map_dfr(all_syns_cart, function(x) {
#   x %$% lm(wgt ~ age + hgt) %$% analyses %>%
#     map_dfr(broom::tidy) %>%
#     group_by(term) %>%
#     summarise(est  = mean(estimate),
#               m    = length(estimate),
#               bm   = sum((estimate - mean(estimate))^2) / (m - 1),
#               ubar = mean(std.error^2),
#               var  = (1 + 1/m) * bm - ubar,
#               pvar = ifelse(var > 0, var, ubar),
#               se   = sqrt(pvar),
#               df   = max(1, (m - 1) * (1 - ubar / ((1 + 1/m) * bm))^2),
#               t    = est / se,
#               p    = pt(t, df))}) %>%
#   group_by(term) %>%
#   summarise(est = mean(est),
#             m   = mean(m),
#             bm = mean(bm),
#             ubar = mean(ubar),
#             var = mean(pvar),
#             se  = mean(se),
#             df  = mean(df),
#             t   = mean(t))

# 
# get_est <- function(mids, formula) {
#   mids %$%
#     lm(wgt ~ age + hgt) %$%
#     analyses %>%
#     map_dfr(broom::tidy) %>%
#     group_by(term) %>%
#     summarise(term = unique(term),
#               est  = mean(estimate),
#               m    = length(estimate),
#               bm   = sum((estimate - mean(estimate))^2) / (m - 1),
#               ubar = mean(std.error^2),
#               var  = (1 + 1/m) * bm - ubar,
#               pvar = ifelse(var > 0, var, ubar),
#               se   = sqrt(pvar),
#               df   = max(1, (m - 1) * (1 - ubar / ((1 + 1/m) * bm))^2),
#               t    = est / se,
#               p    = pt(t, df), .groups = 'drop')
# }
# 
# get_ci <- function(lm_est) {
#   term   <- lm_est$term
#   lower <- lm_est$est + qt(.025, lm_est$df) * lm_est$se
#   upper <- lm_est$est + qt(.975, lm_est$df) * lm_est$se
#   return(tibble(term, lower, upper))
# }
# 
# all_syns_def %>% map_dfr(function(x) {
#   get_est(x, lm(wgt ~ age + hgt)) %>%
#     get_ci %>%
#     mutate(covered = lower < coefs & coefs < upper)
# }) %>% group_by(term) %>% summarise(coverage = mean(covered))
# 
# 
# mice(truth, m = 5, where = random_missings)