# IPM -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# Specify the model
model_file <- "models//ipm//tests//08july2024_testA.txt"
save_file <- "results//ipm//tests//25july_rsltA1.rds"

# Loop dimension parameters
n_year <- 34

# JAGS control parameters
n_i <- 500000
n_a <- 100000
n_b <- 0
n_c <- 1
n_t <- 1000

#Environment====================================================================

require(tidyverse); require(rjags); require(mcmcplots)
# load("data//elk_ipm_data.Rdata")
load("data//elk_ipm_data_14jun2024.Rdata")

#Data_prep======================================================================

p_cjs <- ipm_data$p_cjs
p_cjs[,3] <- p_cjs[,3] + 1
p_cjs <- rbind(p_cjs, ipm_data$p_cjs[1:33,])
p_cjs[67:nrow(p_cjs), 6] <- colSums(ipm_data$calf)[-1]

jags_data <- list(
  s_cjs      = ipm_data$s_cjs,
  p_cjs      = p_cjs,
  # r_ratio    = ipm_data$r_ratio,
  n_sight_ca = ipm_data$n_sight_ca,
  n_sight_am = ipm_data$n_sight_am,
  n_sight_af = ipm_data$n_sight_af,
  n_a_mov    = ipm_data$n_ad_add - abs(ipm_data$n_ad_rem),
  n_c_mov    = ipm_data$n_ca_add - abs(ipm_data$n_ca_rem),
  n_year     = n_year,
  nn_ca      = nrow(ipm_data$n_sight_ca),
  nn_af      = nrow(ipm_data$n_sight_af),
  nn_am      = nrow(ipm_data$n_sight_am),
  nn_count   = nrow(ipm_data$p_cjs),
  ns         = nrow(ipm_data$s_cjs),
  # nr         = nrow(ipm_data$r_ratio),
  n_har      = ipm_data$n_hnt,
  min_ad     = ipm_data$min_ad,
  # min_ca     = ipm_data$min_ca,
  # af_count   = ipm_data$n_f_p_count,
  # nn_fc      = nrow(ipm_data$n_f_p_count),
  # am_count   = ipm_data$n_m_p_count,
  # nn_mc      = nrow(ipm_data$n_m_p_count),
  cdens      = ipm_data$puma_derived,
  n_adj      = ipm_data$af_density,
  min_n1     = ipm_data$min_n1,
  est_n1     = ipm_data$est_n1,
  clim       = ipm_data$palmer_index,
  # y          = ipm_data$y,
  # z          = ipm_data$z,
  # f          = ipm_data$f,
  # l          = ipm_data$l,
  # calf       = ipm_data$calf,
  # male       = ipm_data$male,
  # female     = ipm_data$female,
  # herd       = ipm_data$herd,
  # n_ind      = ipm_data$n_ind,
  n_R_data   = nrow(ipm_data$r_data),
  r_data     = ipm_data$r_data
  # NF_ct      = ipm_data$cjs_n_female,
  # NM_ct      = ipm_data$cjs_n_male
)



inits <- function(){
  # Starting abundance
  N <- array(data = NA, dim = c(4,2,34))
  ca_cjs <- ipm_data$calf * ipm_data$y
  ca_min <- ipm_data$min_ca
  fe_cjs <- ipm_data$female * (1 - ipm_data$calf) * ipm_data$y
  fe_min <- ipm_data$min_ad[1,]
  ma_cjs <- ipm_data$male * (1 - ipm_data$calf) * ipm_data$y
  ma_min <- ipm_data$min_ad[2,]
  for(i in 1:dim(N)[3]){
    N[1,,i] <- (max(ca_cjs[i], ca_min[i]) + 10) / 2
    N[2:4,1,i] <- (max(fe_cjs[i], fe_min[i]) + 10) / 3
    N[2:4,2,i] <- (max(fe_cjs[i], fe_min[i]) + 10) / 3
  }
  N[,,2:34] <- NA
  
  # Recruitment
  R_B0 <- 2
  R_yr <- c(NA, -0.84, -0.55, -1.12, 0.03, -0.68, -0.85, -0.31, 0.00, -0.13, 
            -0.40, -0.59, -0.31, -0.18, -0.26, -0.02, -0.38, -0.99, -0.88, 
            -0.83, -0.49, -0.24, -0.61, -0.78, -0.25, -0.51, -0.36, -0.20,
            -0.33, -0.93, -0.57, -0.52, -0.62, -0.50)
  
  # Survival
  S_C_B0 <- 3
  S_F_B0 <- 3
  S_M_B0 <- 3
  S_C_yr <- c(NA,   0.24,-1.54, 1.61, 1.78, 1.25, 2.11, 1.06, 2.31, 2.16, 2.60, 
              1.71, 2.71, 1.78, 1.63, 1.72, 1.92, 2.32, 0.86, 0.74, 1.14, 1.41, 
              0.78, 0.59, 1.29, 1.13, 1.20, 0.23, 1.73, 1.42, 1.03, 2.13, 0.77, 
              1.38)
  S_F_yr <- c(NA,   1.63, 2.32, 1.36, 3.37, 2.03, 3.15, 2.41, 3.78, 3.08, 2.67, 
              1.82, 3.41, 1.57, 2.15, 2.77, 2.82, 3.18, 1.97, 1.99, 2.87, 2.87,
              2.27, 2.39, 3.22, 2.54, 2.21, 1.98, 3.20, 0.54, 2.59, 2.74, 2.34,
              2.50)
  S_M_yr <- c(NA,   1.39, 0.53, 1.51, 0.42, 1.78, 2.51, 0.68, 1.27, 2.02, 2.07, 
              1.63, 1.90, 1.34, 0.20, 1.65, 1.46, 0.56, 1.12, 1.31, 1.81, 2.03,
              2.42, 1.48, 1.33, 1.40, 1.76, 0.35, 1.02, 1.41, 0.99, 2.43, 1.07, 
              1.45)
  
  # Return
  out <- list(
    N = round(N),
    R_B0 = R_B0,
    S_C_B0 = S_C_B0,
    S_F_B0 = S_F_B0,
    S_M_B0 = S_M_B0,
    R_yr = R_yr,
    S_C_yr = S_C_yr,
    S_M_yr = S_M_yr,
    S_F_yr = S_F_yr
  )
  return(out)
}

initial_values <- inits()

params = c(
  "N_tot",
  "R",
  "S_f",
  "S_m",
  "S_c",
  "N_f",
  "N_c",
  "N_m",
  # "N_lam",
  "R_B0",
  "R_wt",
  "R_wm",
  "R_dd",
  "R_cg",
  "S_C_0",
  "S_wt",
  "S_wm",
  "S_dd",
  "S_cg",
  "lambda",
  "LAMBDA_mean",
  # "pop_r",
  "mean_pop_r",
  "R_mean",
  "S_F_mean",
  "S_M_mean",
  "S_C_mean",
  # "sd_afcount",
  # "sd_amcount",
  "sd_R",
  "S_C_B0",
  "S_F_B0",
  "S_M_B0",
  "sd_S_C",
  "sd_S_M",
  "sd_S_F"
  # "PF",
  # "PM",
  # "P",
  # "exp_ct",
  # "tau_ct",
  # "Ne",
  # "Te",
  # "Naug"
)

#Model==========================================================================
cat("This run was started at", as.character(Sys.time()), "\n")
tictoc::tic()

jgs_mdl <- jags.model(
  file = model_file,
  data = jags_data,
  inits = inits,
  n.chains = n_c,
  n.adapt = n_a
)

# update(jgs_mdl, n.iter = n_i*2)

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

