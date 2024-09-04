# Run the cjs model
# Kenneth Loonam
# August 2024

# Load packages
require(nimble); require(lubridate)

# Specify results location
result_file <- paste0(
  "results//cjs_", 
  today(), 
  ".rds"
)

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
  "prob_am"
)

# Run it
rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = params,
  inits       = cd$initial_values,
  niter       = 250,
  nburnin     = 100,
  nchains     = 3,
  progressBar = T,
  summary     = F,
  check       = F
)  

save(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
