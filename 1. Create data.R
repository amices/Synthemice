# Create data sets to work with

# Boys data

# Stappenplan:
# 1. partitioneer de boys data in 5 delen (748/5 +/- 150 per data set)
# 2. maak synthetic versions van de 5 delen
# 3. Merge de 5 synthetic delen in een centrale synthetic data set
# 4. Imputeer de centrale synthetic data set
# 5. Deconstruct de m imputed sets into the 5 synthetic sets
# 6. Match the corresponding synthetic imputation to the observed imputation in 
# all of the (m by 5) synthetic imputed sets

# Wanneer zijn we tevreden?
# - Has algorithm converged (obv synthetic version --> evt later obv de real data version)
# - Hanne's suggestion for tracking a model/data parameter based on the real data
# - Als bias, CI width en Coverage of the 95% CI about the parameter (e.g. beta or mean) ok is. 

# Hoe gaan we te werk
#
# - Proof of Concept: Incomplete boys data als vertrekpunt en kijken:
#       1. is geconvergeerd?
#       2. is regressiemodel obv federated anders dan regressiemodel obv unfederated. 
#       3. 1000 sims. 
#
# - Als PoC succesvol --> jeej
#       1. Full blown sim met complete(mice(boys, m = 1, seed = 123)) als TRUE/population input
#       2. Bias? CI width? Coverage?

## Packages
require(mice)
require(tidyverse)
require(modelr)

## 
set.seed(123)

## Partition data
ngroups <- 5                                    # specify number of subsets
partition <- rep(1/ngroups,ngroups)             # create partitioning object
names(partition) <- paste0("data", 1:ngroups)   # add names (data1 to data5)
part <- resample_partition(boys, partition)     # create the partitioned data

#part %>%                                        # use map function to perform lm on all 5 subsets
#  map(., ~ lm(bmi ~ age + hgt, data = .))       # just to illustrate the use of map for myself


