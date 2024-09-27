# IPM -- Starkey Elk
# Kenneth Loonam
#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)

#Data===========================================================================
# Prepare the full IPM data
year_range <- 1988:2023
max_age <- 2
climate_cov <- "spei_3m" 
  # Choose from:
    # pdi_growing, pdi_september, temp, precip, ndvi, or spei_(1,3,6,12,24,48)m
puma_cov <- "pd_full_mean" 
  # Choose from:
    # pd_reconstruction, pd_mortalities, pd_odfwe_est, pd_logistic, pd_full_mean
    # and pd_trim_mean
source("workflows//data_prep//ipm_data_bundling.R")
# Prepare the IPM initial values
source("workflows//data_prep//ipm_build_inits.R")

#Variables======================================================================
# Specify the model and save location for the results
model_file <- "models//ipm_1.1.txt"
save_file <- "results//submission_2//ipm_1.0_fullmean_spei3m.rds"

# JAGS control parameters
n_i <- 1000000
n_a <- 10000
n_b <- 100000
n_c <- 3
n_t <- 100

#Data_prep======================================================================
jags_data <- readRDS("data//ipm_data.rds")
jags_inits <- readRDS("data//ipm_inits.rds")
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
