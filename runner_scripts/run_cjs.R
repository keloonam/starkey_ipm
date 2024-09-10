# Run the cjs model
# Kenneth Loonam
# August 2024

# Load packages
require(nimble); require(lubridate)

# Specify results location
result_file <- "results//cjs_rslt_06sep2024.rds"

# Load the data
cd <- readRDS("data//cjs_data.rds")

# Load the model
source("models//cjs.R")

# Set trace monitors
params <- c(
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "prob_af",
  "prob_am",
  "post_diff_mn",
  "post_diff_sd"
)

# Run it
rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = params,
  inits       = cd$initial_values,
  niter       = 25000,
  nburnin     = 10000,
  nchains     = 3,
  progressBar = T,
  summary     = F,
  check       = F
)  

saveRDS(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
