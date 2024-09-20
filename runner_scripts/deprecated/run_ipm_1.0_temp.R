# IPM -- Starkey Elk
# Kenneth Loonam
#Variables======================================================================

# Specify the model and save location for the results
model_file <- "models//ipm_1.0.txt"
save_file <- "results//ipm_1.0_rslts_17sep2024.Rdata"

# JAGS control parameters
n_i <- 100000
n_a <- 1000
n_b <- 10000
n_c <- 3
n_t <- 10

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)

#Data_prep======================================================================

jags_data <- readRDS("data//ipm_data.rds")
jags_data$cjs_c <- jags_data$cjs_c %>% as_tibble() %>%
  mutate(cnt = count / prob) %>%
  select(class, yr, cnt) %>%
  as.matrix()
jags_inits <- readRDS("data//ipm_inits.rds")
params = c(
  "N", "NF", "NM", "NC", "Ntot",
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

# update(jgs_mdl, n.iter = n_b)

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

res <- rslt %>% map(as_tibble) %>%
  bind_rows() 
res %>%
  select(grep("NC", names(.)))
res %>% select(grep("Ntot", names(.)))
