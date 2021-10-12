# Estimate elk harvest rate by year/age/sex
# Kenneth Loonam
# October 2020

# On second thought, I'm not sure this is independently estimable
# This script currently does nothing. 
# Resuming work on the survival-harvest approach.

#Variables======================================================================

start_year <- 1988
end_year <- 2020

# Parameters to track
params <- c(

)

# Tau for logistic transformed priors
pr_p <- 0.5

# Model
model_file <- "models//harvest//harvest_model.txt"

# Sampler variables
ni <- 100
nt <- 1
nb <- 20
nc <- 3

#Environment====================================================================

require(tidyverse); require(nimble)
load("data//elk_data.Rdata")

#Data_prep======================================================================




