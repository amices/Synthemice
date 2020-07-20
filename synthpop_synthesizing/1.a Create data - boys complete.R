###############################################################################
## Create data sets to work with, we will work with the boys dataset that is ##
## build-in in mice.                                                         ##
###############################################################################


###############################################################################
## Load packages and set seed                                                ##
###############################################################################

library(mice)   ## mice for the boys data and the imputations
library(modelr) ## modelr to do the partitioning   

set.seed(123)   # seed for reproducibility of the imputations and partitioning

###############################################################################
## obtain and impute the boys dataset, create the completed boys data, and   ##
## partition the completed dataset                                           ##
###############################################################################

meth <- make.method(boys)
pred <- make.predictorMatrix(boys)
meth["bmi"] <- "~ I(wgt / (hgt/100)^2)"

pred[c("wgt", "hgt"), "bmi"] <- 0

imp <- mice(boys, meth = meth, pred = pred,      ## Impute the boys dataset by means of a single
            m = 1, printFlag = F)                ## imputation, and fill in the imputed data into
boyscomp <- complete(imp)                        ## the dataset containing missing values.
                     

ngroups <- 5                                     ## Specify five groups, and
partitions <- rep(1/ngroups,ngroups)             ## specify five equal fractions
names(partitions) <- paste0("d", 1:ngroups)      ## give names to the partitionings and
part <- resample_partition(boyscomp, partitions) ## partition the data in five equal parts
