# Run the cjs model
# Kenneth Loonam
# August 2024

# Load packages
require(nimble); require(lubridate)

# Specify results location
result_file <- "results//recruitment_rslt_06sep2024.rds"

# Load the data
cd <- readRDS("data//recruitment_data.rds")

# Load the model
source("models//recruitment.R")

# Set trace monitors
params <- c(

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
