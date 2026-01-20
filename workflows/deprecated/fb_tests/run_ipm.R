# IPM -- Starkey Elk
# Kenneth Loonam
#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)

#Data===========================================================================
# Prepare the full IPM data
yr_range <- 2018:2023
source("fb_tests//ipm_data_bundling.R")
# Prepare the IPM initial values
source("fb_tests//ipm_build_inits.R")

#Variables======================================================================
# Specify the model and save location for the results
model_file <- "fb_tests//ipm_fb.txt"
yr_name <- paste0(yr_range[1], "_to_", yr_range[length(yr_range)])
save_file <- paste0("fb_tests//results//fb_", yr_name, "_rslt.rds")

# JAGS control parameters
n_i <- 500000
n_a <- 10000
n_b <- 250000
n_c <- 3
n_t <- 50

#Data_prep======================================================================
jags_data <- readRDS("fb_tests//ipm_data.rds")
jags_inits <- readRDS("fb_tests//ipm_inits.rds")
params = c(
  "NF", "NM", "NC", "Ntot", "sd_c",
  "LAMBDA", "MEAN_LAMBDA",
  "R", "R_B0", "R_WT", "R_WM", "R_CG", "R_DD", "sd_R",
  "SC", "SC_B0", "SC_WT", "SC_WM", "SC_CG", "SC_DD", "sd_SC",
  "SF", "SM", "sd_SF", "sd_SM"
)

#Model==========================================================================
cat("This run was started at", as.character(Sys.time()), "\n")
tictoc::tic()

jgs_mdl <- jags.model(
  file = model_file,
  data = jags_data,
  inits = jags_inits,
  n.chains = n_c,
  n.adapt = n_a
)

update(jgs_mdl, n.iter = n_b)

rslt <- coda.samples(
  jgs_mdl,
  variable.names = params,
  n.iter = n_i,
  thin = n_t
)

cat("From start to end of this model run,")
tictoc::toc()
mcmcplots::mcmcplot(rslt)
saveRDS(rslt, file = save_file)
