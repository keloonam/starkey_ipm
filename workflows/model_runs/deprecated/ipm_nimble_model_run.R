# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//ipm_elk_null_boar_style.R"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 6000
n_b <- 1000
n_c <- 3
n_t <- 10

#Environment====================================================================

require(tidyverse); require(mcmcplots); require(nimble)
load("data//elk_ipm_nimble_data.Rdata")

