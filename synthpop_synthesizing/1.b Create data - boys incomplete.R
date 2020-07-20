###############################################################################
## Create data sets to work with, we will work with the boys dataset that is ##
## build-in in mice.                                                         ##
###############################################################################


###############################################################################
## Load packages and set seed                                                ##
###############################################################################

require(mice)   ## mice for the boys data and the imputations
require(modelr) ## modelr to do the partitioning   

set.seed(123)   # seed for reproducibility of the imputations and partitioning

###############################################################################
## Obtain the boys dataset, and partition it into five equal parts           ##
###############################################################################

## Partition data
ngroups <- 5                                    # specify number of subsets
partition <- rep(1/ngroups,ngroups)             # create partitioning object
names(partition) <- paste0("data", 1:ngroups)   # add names (data1 to data5)
part <- resample_partition(boys, partition)     # create the partitioned data

pool(part %>%                                        # use map function to perform lm on all 5 subsets
       map(., ~ lm(bmi ~ age + hgt, data = .)))       # just to illustrate the use of map for myself

lm(bmi ~ age + hgt, data = boys)

