# IPM -- Starkey Elk
# Kenneth Loonam
#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)

#Data===========================================================================
# Prepare the full IPM data
# yr_range <- 1988:1994
# yr_range <- 1994:2000
# yr_range <- 2000:2006
# yr_range <- 2006:2012
# yr_range <- 2012:2018
yr_range <- 2018:2023
source("fb_tests//ipm_data_bundling.R")
# Prepare the IPM initial values
source("fb_tests//ipm_build_inits.R")

#Variables======================================================================
# Specify the model and save location for the results
model_file_fb <- "fb_tests//ipm_fb.txt"
yr_name <- paste0(yr_range[1], "_to_", yr_range[length(yr_range)])
save_file_fb <- paste0("fb_tests//results//fb_", yr_name, "_rslt.rds")

# JAGS control parameters
n_i_fb <- 500000
n_a_fb <- 10000
n_b_fb <- 25000
n_c_fb <- 3
n_t_fb <- 50

n_i_pb <- 500000
n_a_pb <- 10000
n_b_pb <- 100000
n_c_pb <- 3
n_t_pb <- 50

cjs_ni <- 50000
cjs_nb <- 10000

params = c(
  "NF", "NM", "NC", "Ntot", "sd_c",
  "LAMBDA", 
  "R", "R_B0", "R_WT", "R_WM", "R_CG", "R_DD", "sd_R",
  "SC", "SC_B0", "SC_WT", "SC_WM", "SC_CG", "SC_DD", "sd_SC",
  "SF", "SM", "sd_SF", "sd_SM"
)

#Model==========================================================================
cat("This run was started at", as.character(Sys.time()), "\n")
tictoc::tic()

jgs_mdl_fb <- jags.model(
  file = model_file_fb,
  data = jags_data,
  inits = jags_inits,
  n.chains = n_c_fb,
  n.adapt = n_a_fb
)

update(jgs_mdl_fb, n.iter = n_b_fb)

rslt_fb <- coda.samples(
  jgs_mdl_fb,
  variable.names = params,
  n.iter = n_i_fb,
  thin = n_t_fb
)

cat("From start to end of this model run,")
tictoc::toc()
saveRDS(rslt_fb, file = save_file_fb)

require(nimble)
source("models//survival.R")

# Set trace monitors
cjs_params <- c(
  "survival_af", 
  "survival_am", 
  "survival_ca",
  "prob_af",
  "prob_am",
  "pp_mean_difference",
  "pp_sd_difference"
)

cd <- list(
  constants = list(
    nocc = ncol(jags_data$y),
    nind = nrow(jags_data$y),
    l = jags_data$l,
    f = jags_data$f
  ),
  data = list(
    y = jags_data$y,
    z = jags_data$z,
    m = jags_data$m,
    h = jags_data$h,
    c = jags_data$c
  ),
  inits = list(
    bsm = rep(0, ncol(jags_data$y)),
    bsc = rep(0, ncol(jags_data$y)),
    bpm = rep(0, ncol(jags_data$y)),
    bsh = rep(0, ncol(jags_data$y)),
    bph = rep(0, ncol(jags_data$y))
  )
)

# Run it
cjs_rslt <- nimbleMCMC(
  code        = code,
  constants   = cd$constants,
  data        = cd$data,
  monitors    = cjs_params,
  inits       = cd$inits,
  niter       = cjs_ni,
  nburnin     = cjs_nb,
  nchains     = 3,
  thin        = 5,
  progressBar = T,
  summary     = F,
  check       = F
)  

logit <- function(x){
  log(x/(1-x)) %>% return()
}

rm_bad_yrs <- function(x){
  x %>%
    filter(yr != "2017") %>%
    filter(yr != "2020") %>%
    return()
}

cjs_rslt <- cjs_rslt %>%
  map(as_tibble) %>%
  bind_rows()

#Capture Probability============================================================
fp_logit <- cjs_rslt %>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bpf")
mp_logit <- cjs_rslt %>%
  select(grep("prob_am", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bpm")

#Survival=======================================================================
fs_logit <- cjs_rslt %>%
  select(grep("survival_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsf")
ms_logit <- cjs_rslt %>%
  select(grep("survival_am", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsm")
cs_logit <- cjs_rslt %>%
  select(grep("survival_ca", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsc")

#Cleanup Steps==================================================================
survival_summary <- list(
  sf_logit = fs_logit,
  sm_logit = ms_logit,
  sc_logit = cs_logit,
  pf_logit = fp_logit,
  pm_logit = mp_logit
)


survival_summary <- map(survival_summary, rm_bad_yrs)

est_bs <- bind_rows(survival_summary) %>%
  filter(yr %in% yr_range[-c(1, length(yr_range)-1, length(yr_range))]) %>%
  filter(var %in% c("bsf", "bsc", "bsm")) %>%
  mutate(yr = as.numeric(yr) - (yr_range[1] - 1)) %>%
  mutate(class = case_when(
    var == "bsf" ~ 2,
    var == "bsc" ~ 1,
    var == "bsm" ~ 3
  )) %>%
  mutate(tau = sd_to_tau(sd)) %>%
  select(class, yr, mn, tau) %>%
  arrange(yr, class) %>%
  as.matrix()

est_br <- readRDS("results//recruitment_summary.rds") %>%
  filter(yr %in% yr_range) %>%
  mutate(class = 0) %>%
  mutate(yr = as.numeric(yr) - (yr_range[1] - 1)) %>%
  mutate(tau = sd_to_tau(sd)) %>%
  select(class, yr, mn, tau) %>%
  arrange(yr) %>%
  filter(yr > 1) %>%
  as.matrix()

jags_data$est_br <- est_br
jags_data$est_bs <- est_bs
jags_data$ns <- nrow(est_bs)
jags_data$nr <- nrow(est_br)

save_file_pb <- paste0("fb_tests//results//pb_", yr_name, "_rslt.rds")
#Data_prep======================================================================

#Model==========================================================================
cat("This run was started at", as.character(Sys.time()), "\n")
tictoc::tic()

jgs_mdl_pb <- jags.model(
  file = "models//ipm_1.0.txt",
  data = jags_data,
  inits = jags_inits,
  n.chains = n_c_pb,
  n.adapt = n_a_pb
)

update(jgs_mdl_pb, n.iter = n_b_pb)

rslt_pb <- coda.samples(
  jgs_mdl_pb,
  variable.names = params,
  n.iter = n_i_pb,
  thin = n_t_pb
)

cat("From start to end of this model run,")
tictoc::toc()
saveRDS(rslt_pb, file = save_file_pb)
