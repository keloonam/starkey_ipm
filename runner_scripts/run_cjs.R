# Run the cjs model
# Kenneth Loonam
# August 2024

# Load packages
require(nimble); require(lubridate)

# Specify results location
result_file <- "results//cjs_rslt.rds"

# Load the data
cd <- readRDS("data//cjs_data.rds")

# Load the model
source("models//survival.R")

# Set trace monitors
params <- c(
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "prob_af",
  "prob_am",
  "pp_mean_difference",
  "pp_sd_difference"
)

# Run it
rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = params,
  inits       = cd$initial_values,
  niter       = 50000,
  nburnin     = 15000,
  nchains     = 3,
  thin        = 5,
  progressBar = T,
  summary     = F,
  check       = F
)  

saveRDS(rslt, file = result_file)
# mcmcplots::mcmcplot(rslt)
rm(list = ls()[-which(ls() == "yr_range")])
