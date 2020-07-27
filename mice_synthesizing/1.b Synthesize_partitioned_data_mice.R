


n.parts <- 5

parts <- rep(1/n.parts, n.parts)
names(parts) <- paste0("P", 1:n.parts)

boot_part <- bootstrap_boys %>%
  future_map(function(x) resample_partition(x, parts) %>%
               map(function(y) as.data.frame(y)))