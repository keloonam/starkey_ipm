# IPM run - Starkey Elk
# Kenneth Loonam
# May 2024
results_file_name <- "results//ipm_24jun2024_test.rds"
#Environment====================================================================

require(nimble); require(mcmcplots); require(dplyr); require(rjags)

full_data <- readRDS("data//the_ipm_data.rds")

source("models//ipm//the_ipm.R")
load("data//elk_ipm_data_21apr2023.Rdata")

nchains <- 3
nburnin <- 10000
niter <- 30000
nthin <- 5
#Data===========================================================================

dtf <- list(
  NC_est = full_data$NC_est[,2:3],
  NF_est = full_data$NF_est[,2:3],
  NM_est = full_data$NM_est[,2:3],
  r_dt = full_data$r_dt[,2:3],
  nc1e = full_data$nc1e,
  nf1e = full_data$nf1e,
  nm1e = full_data$nm1e,
  # SCcons = rep(1, 34),
  # SFcons = rep(1, 34),
  # SMcons = rep(1, 34),
  # Rcons  = rep(1, 34),
  # AugFmin = abs(full_data$NFmin),
  # AugMmin = abs(full_data$NMmin),
  # NCmin = full_data$NCmin,
  # y = full_data$y,
  # z = full_data$z,
  # NF_ct = full_data$NF_ct,
  # NM_ct = full_data$NM_ct
  af_count = ipm_data$n_f_p_count[,4],
  am_count = ipm_data$n_m_p_count[,4],
  sc_cjs = ipm_data$s_cjs[62:92,4:5],
  sf_cjs = ipm_data$s_cjs[1:30, 4:5],
  sm_cjs = ipm_data$s_cjs[31:61,4:5]
)

#Constants======================================================================
cnst <- list(
  NFman = full_data$NFman,
  NCman = full_data$NCman,
  NMman = full_data$NMman,
  NFhar = full_data$NFhar,
  NMhar = full_data$NMhar,
  n_year = full_data$n_year,
  nNC = full_data$nNC,
  nNF = full_data$nNF,
  nNM = full_data$nNM,
  nR = full_data$nR,
  nc1e_min = full_data$nc1e_min,
  nf1e_min = full_data$nf1e_min,
  nm1e_min = full_data$nm1e_min,
  veg = full_data$sep_pdi,
  pum = full_data$puma_composit,
  elk = full_data$elk_density,
  NC_est_t = full_data$NC_est[,1],
  NF_est_t = full_data$NF_est[,1],
  NM_est_t = full_data$NM_est[,1],
  r_dt_t = full_data$r_dt[,1],
  # l = full_data$l,
  # male = full_data$male,
  # female = full_data$female,
  # calf = full_data$calf,
  # herd = full_data$herd,
  # f = full_data$f,
  # nind = full_data$nind
  nn_fc = length(dtf$af_count),
  nn_mc = length(dtf$am_count),
  af_count_t = ipm_data$n_f_p_count[,1],
  am_count_t = ipm_data$n_f_p_count[,1],
  nsc = nrow(dtf$sc_cjs),
  nsf = nrow(dtf$sf_cjs),
  nsm = nrow(dtf$sm_cjs),
  sc_cjs_t = ipm_data$s_cjs[62:92,1],
  sf_cjs_t = ipm_data$s_cjs[1:30, 1],
  sm_cjs_t = ipm_data$s_cjs[31:61,1]
)

#Initial values=================================================================

inits <- readRDS("data//the_ipm_inits.rds")
inits <- inits[c(1:31, 40:46)]
inits$sd_amcount <- 25
inits$sd_afcount <- 25

#Model setup====================================================================
mons <- c(
  "R", "SC", "SF", "SM", 
  "NC", "NF", "NM", 
  "R_Bvy", "R_Bvm", "R_B0_mean", "R_Bpu", "R_Bdd", "R_B0_sd", 
  "SC_B0_mean", "SC_B0_sd", "SC_Bvy", "SC_Bvm", "SC_Bpu", "SC_Bdd", 
  "NFaug", "NMaug", "NCaug", 
  "SM_B0_mean", "SF_B0_mean", "SM_B0_sd", "SF_B0_sd", 
  "Ntot", "LAMBDA"
)
# "PF", "PM",

#Run the IPM in nimble==========================================================
# # This doesn't work (yet?)
# # Look, models like this are hard, and nimble is complicated
# # I really want to revisit this when I have infinite time


ipm <- nimbleModel(
  code = nimble_code,
  constants = cnst,
  data = dtf,
  inits = inits
)

CMipm <- compileNimble(ipm)

cnf <- configureMCMC(
  model = ipm,
  monitors = mons,
  onlySlice = F,
  autoBlock = F
)

MCipm <- buildMCMC(ipm)
MCMCipm <- compileNimble(MCipm)

rslt <- runMCMC(
  MCMCipm,
  niter = 50000,
  nburnin = 0,
  nchains = 3,
  inits = inits,
  thin = 10
)

mcmcplot(rslt)
