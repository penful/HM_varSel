#########################################################
#### Code provided for the scientific paper titled: 
### "Variable selection for hidden Markov models with continuous variables and missing data" 
### Authors: F. Pennoni, F. Bartolucci, F. Pandolfi 
### Journal of Classification (2023)
#########################################################

#### Example to run #####
### It requires almost 2 minutes  of execution time ###
rm(list=ls())
require(LMest)
library(mvtnorm)
require(Formula)
require(MASS)
require(mix)

source("regress_miss.R")         #Function to fit the multivariate regression model
source("lmestContMISS.R")        #Function to estimate the HM model with missing data
source("functions.R")            #Internal function
source("lmbasic.cont.MISS.R")    #Internal function
source("complk_cont_miss.R")     #Internal function
source("drawHMBasicCont.R")
source("compute_BIC.R")
source("item_selection.R")
source("count_eq.R")
source("forward_regress_miss.R")

# load data of countries #
load("dt1.Rdata") 

# specify the option for time heterogeneous transition probabilities
modBasic <- 0 

# print the required time and run the main function to perform variables and states selection 
start_time <- Sys.time()
out1 <-  item_selection(Y = dt1,
                        Kmax = 5,
                        modBasic=modBasic)
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

# Print the execution time in minutes
print(execution_time)
# Save the output file
save.image("example_trial.Rdata")

