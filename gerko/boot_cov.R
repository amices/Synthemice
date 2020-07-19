# parameters
set.seed(123)

# packages
library(mvtnorm)
library(mice)
library(plyr)
library(dplyr)
library(magrittr)
library(purrr)

set.seed(123)
nsim = 10000

# establish origin
truth <- mice(boys, m=1, print = FALSE) %>% mice::complete()
truemodel <- truth %$%
  lm(wgt ~ hgt + age)

# sim bootstrap
simdata <- replicate(nsim, 
                     truth[sample(1:748, 748, replace = TRUE), ], 
                     simplify = FALSE)

# calculate statistics
model <- simdata %>% map(~ lm(wgt ~ hgt + age, data = .x))
ci <- model %>% map(confint)
estimates <- model %>% 
  map(coef) %>% 
  do.call(rbind, args = .)
cov <- ci %>% 
  map(function(x) x[, 1] <= coef(truemodel) & coef(truemodel) <= x[, 2]) %>%
  do.call(rbind, args = .)
manual <- model %>% 
  map(function(x) cbind(vcov(x) %>% diag %>% sqrt(.), coef(x))) %>% # est and resp. vars
  map(function(x) cbind(x[, 2] - qt(.975, 747) * x[, 1],  x[, 2] + qt(.975, 747) * x[, 1])) %>%
  map(function(x) x[, 1] <= coef(truemodel) & coef(truemodel) <= x[, 2]) %>%
  do.call(rbind, args = .)
manual2 <- model %>% 
  map(function(x) cbind(sqrt(diag(var(estimates))), coef(x))) %>% # sd of 
  map(function(x) cbind(x[, 2] - qt(.975, 747) * x[, 1],  x[, 2] + qt(.975, 747) * x[, 1])) %>%
  map(function(x) x[, 1] <= coef(truemodel) & coef(truemodel) <= x[, 2]) %>%
  do.call(rbind, args = .)

# manual based on empirical quantiles
intercept <- data.frame(est = mean(estimates[, 1]), 
                        ciw = diff(quantile(estimates[, 1], probs = c(.025, .975)))) %>%
  mutate(low = est - ciw/2,
         up = est + ciw/2)
hgt <- data.frame(est = mean(estimates[, 2]), 
                  ciw = diff(quantile(estimates[, 2], probs = c(.025, .975)))) %>%
  mutate(low = est - ciw/2,
         up = est + ciw/2)
age <- data.frame(est = mean(estimates[, 3]), 
                  ciw = diff(quantile(estimates[, 3], probs = c(.025, .975)))) %>%
  mutate(low = est - ciw/2,
         up = est + ciw/2)
combined <- rbind(intercept, hgt, age)

covint <- data.frame(est = estimates[, 1], 
                     low = estimates[, 1] - intercept$ciw/2,
                     up = estimates[, 1] + intercept$ciw/2) %>%
  mutate(cov = low <= coef(truemodel)[1] & coef(truemodel)[1] <= up)
covhgt <- data.frame(est = estimates[, 2], 
                     low = estimates[, 2] - hgt$ciw/2,
                     up = estimates[, 2] + hgt$ciw/2) %>%
  mutate(cov = low <= coef(truemodel)[2] & coef(truemodel)[2] <= up)
covage <- data.frame(est = estimates[, 3], 
                     low = estimates[, 3] - age$ciw/2,
                     up = estimates[, 3] + age$ciw/2) %>%
  mutate(cov = low <= coef(truemodel)[3] & coef(truemodel)[3] <= up)


# bias
colMeans(estimates) - coef(truemodel)

# coverage
colMeans(cov) # wrong
colMeans(manual) # identical, but wrong
colMeans(manual2) # correct - empirical
colMeans(cbind(covint$cov, covhgt$cov, covage$cov)) # correct - empirical 

