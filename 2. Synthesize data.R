# Create synthetic versions of the data

# Obtain the previously synthesized data
source("1. Create data.R")

##### This is all very preliminary.

syns <- lapply(part, function(x) syn(x[[1]]))

syns$group1$syn

part %>%
  map(., ~ syn(.))

group1 <- part$group1[[1]]

install.packages("synthpop")

syn1 <- syn(group1)

syn1$syn
summary(syn1$syn)
summary(group1)
compare(syn1, group1)