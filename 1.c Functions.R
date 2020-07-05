#########################################################################################
## Function for synthesizing the data, that handles partitioned data and the           ##
## complete data                                                                       ##
#########################################################################################

synthesize <- function(data, parts = NULL,    # Use a function to synthesize and analyze
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
    syn_dat <- map_dfr(parts, function(x) syn(as.data.frame(x), print.flag = F, m = m, method = method)$syn) 
  }                                           
  else {
    # If it concerns the complete data as a whole, simply synthesize the complete data, 
    # and extract this data
    if (m == 1) syn_dat <- syn(data, print.flag = F, m = m, method = method)$syn
    else syn_dat <- syn(data, print.flag = F, m = m, method = method)
  }
  return(syn_dat)
}

#########################################################################################
## Function for synthesizing the data, that handles partitioned data and the           ##
## complete data                                                                       ##
#########################################################################################

normal_syn <- function(data, method = NULL, formula, pop.inf = F) {
  d   <- data                                             # specify the data
  s1  <- syn(d, m = 1, method = method, print.flag = F)   # 1 synthesis
  s5  <- syn(d, m = 5, method = method, print.flag = F)   # 5 syntheses
  s10 <- syn(d, m = 10, method = method, print.flag = F)  # 10 syntheses
  
  f <- as.formula(formula)                                # the formula
  
  # Fit the lm models on all four datasets (the real data and the synthetic data)
  r <- list(Real_Sample = lm(f,d), Syn1 = lm.synds(f, s1), Syn5 = lm.synds(f, s5), Syn10 = lm.synds(f, s10))
  
  # Collect the output
  out <- map(r, function(x) {if (class(x) == "fit.synds") {x <- summary(x, population.inference = pop.inf)}
    else {x <- summary(x)}
    return(as.data.frame(coef(x)))}) %>% 
    map(., function(x) {colnames(x) <- c("Est", "SE", "stat", "P"); return(as.data.frame(x))}) %>%
    map(., function(x) rownames_to_column(x, var = "Variable")) %>%
    bind_rows(., .id = "Method") %>%
    group_by(Method) %>% 
    ungroup() %>%
    dplyr::select(., Method, Variable, Est, SE, stat, P) %>%
    mutate(Lower = Est - qnorm(.975)*SE,
           Upper = Est + qnorm(.975)*SE)
  # Return the output
  return(as.data.frame(out))
}

# Generate multicariate normal data
normal <- function(r2, ratio_beta, rho, n = 10000) {
  X <- MASS::mvrnorm(n = n, mu = rep(0, length(ratio_beta)), Sigma = rho)
  colnames(X) <- paste0("IV", 1:length(ratio_beta))
  coefs <- matrix(NA, nrow = length(ratio_beta), ncol = length(ratio_beta))
  for (i in 1:length(ratio_beta)) {
    for (j in 1:length(ratio_beta)) {
      coefs[i,j] <- ifelse(i > j, ratio_beta[i] * ratio_beta[j], 0)
    }
  }
  b <- sqrt((r2 / (sum(ratio_beta^2) + 2 * sum(coefs * rho))))
  y <- X %*% (b*ratio_beta) + rnorm(n, 0, sqrt(1 - r2))
  return(list(dat = data.frame(DV = y, X), coefs = b*ratio_beta))
}

syn_part_data <- function(partitioned_data, M, formula, 
                          method = NULL, pop.inf = T) {
  
  d <- partitioned_data
  f <- formula
  
  out <- map(d, function(x) syn(as.data.frame(x), m = M, print.flag = F,
                                method = method))
  
  new_syns <- as.list(1:M)
  
  for (m in 1:M) {
    new_syns[[m]] <- map_dfr(1:length(out), function(x) out[[x]]$syn[[m]])
  }
  
  syn <- out[[1]]
  syn$syn <- new_syns
  
  data <- map_dfr(d, function(x) as.data.frame(x))
  
  # Fit the lm models on all datasets
  r <- list(Real_Sample = lm(f,data), Syn = lm.synds(f, syn))
  
  # Collect the output
  out <- map(r, function(x) {if (class(x) == "fit.synds") {x <- summary(x, population.inference = pop.inf)}
    else {x <- summary(x)}
    return(as.data.frame(coef(x)))}) %>% 
    map(., function(x) {colnames(x) <- c("Est", "SE", "stat", "P"); return(as.data.frame(x))}) %>%
    map(., function(x) rownames_to_column(x, var = "Variable")) %>%
    bind_rows(., .id = "Method") %>%
    group_by(Method) %>% 
    ungroup() %>%
    dplyr::select(., Method, Variable, Est, SE, stat, P) %>%
    mutate(Lower = Est - qnorm(.975)*SE,
           Upper = Est + qnorm(.975)*SE)
  # Return the output
  return(as.data.frame(out))
}

syn_part_data_new_var <- function(partitioned_data, M, formula, pop.inf = T) {
  
  d <- partitioned_data
  f <- formula
  
  out <- map(d, function(x) syn(as.data.frame(x), m = M, print.flag = F))
  
  all_syns <- as.list(1:M)
  
  for (m in 1:M) {
    all_syns[[m]] <- map(1:length(out), function(x) out[[x]]$syn[[m]])
  }
  
  results <- all_syns %>% unlist(recursive = F) %>% 
    map_dfr(., .id = "syn", function(x) {
      lm(f, x) %>% broom::tidy() %>%
      mutate(variance = std.error^2) %>%
      dplyr::select(term = term, est = estimate, variance = variance)})
  
  real_results <- map_dfr(d, function(x) as.data.frame(x)) %>%
                  lm(f, .) %>% broom::tidy(conf.int=T) %>%
                  dplyr::select(term = term, Qbar = estimate, SE = std.error,
                                Lower = conf.low, Upper = conf.high)
  
  results %>%
    group_by(term) %>%
    summarise(Qbar = mean(est),
              bm = sum((est - Qbar)^2 / (n()-1)),
              vbar = mean(variance),
              Ts = (bm * (1 + 1 / (n())) -vbar),
              Tstar_s = max(0,Ts) + ifelse(Ts < 0, 1, 0)*vbar,
              SE = sqrt(Tstar_s),
              rm = (1 + 1/n())*bm/vbar,
              df = (n() - 1)*(1 - 1/rm)^2,
              Lower = Qbar - qt(.975, df)*sqrt(Tstar_s),
              Upper = Qbar + qt(.975, df)*sqrt(Tstar_s)) %>% ungroup() -> summary_results
    #dplyr::select(term, Qbar, SE, Lower, Upper) -> summary_results
  
  output <- bind_rows(Real = real_results, Synthetic = summary_results, .id = "Method")
  
  return(output)
}

print_results <- function(output, coefs) {
  output %>%
    mutate(RealEst = rep(coefs, nrow(output)/length(coefs)),
           Covered = Lower < coefs & Upper > coefs,
           Bias = Est - RealEst) %>%
    ungroup() %>% group_by(Method, Variable) %>%
    summarise("Population estimate" = unique(RealEst),
              "Qbar" = mean(Est),
              "Bias" = mean(Bias),
              "MeanSE" = mean(SE),
              "MinSE" = min(SE),
              "MaxSE" = max(SE),
              "Lower" = mean(Lower),
              "Upper" = mean(Upper),
              "CIW" = mean(Upper) - mean(Lower),
              "Coverage" = mean(Covered)) %>%
    arrange(factor(Method, levels = c("Real_Sample", "Syn1", "Syn5", "Syn10")))
}

