# Run the cjs model
# Kenneth Loonam
# August 2024

# Load packages
require(nimble); require(lubridate)

# Specify results location
result_file <- "results//mscr_rslt.rds"

# Load the data
cd <- readRDS("data//cjs_data.rds")

# Load the model
source("models//mscr.R")

# Set trace monitors
params <- c(
  "Pf", "Pm", 
  "Snaf", "Snam", "Snca",
  "Shaf", "Sham", "Shjf", "Shjm",
  "Saf", "Sam", "Sjf", "Sjm",
  "fit.new", "fit.org"
)

# Run it
rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = params,
  inits       = cd$initial_values,
  niter       = 20000,
  nburnin     = 5000,
  nchains     = 3,
  thin        = 3,
  progressBar = T,
  summary     = F,
  check       = F
)  

saveRDS(rslt, file = result_file)
mcmcplots::mcmcplot(rslt)
rm(list = ls()[-which(ls() == "yr_range")])
