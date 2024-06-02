# IPM run - Starkey Elk
# Kenneth Loonam
# May 2024

#Environment====================================================================

require(nimble); require(mcmcplots); require(dplyr)

full_data <- readRDS("data//the_ipm_data.rds")

source("models//ipm//the_ipm_less_cjs.R")
load("data//elk_ipm_data_21apr2023.Rdata")

#Data===========================================================================

dtf <- list(
  # y = full_data$y,
  # z = full_data$z,
  NC_est = full_data$NC_est[,2:3],
  NF_est = full_data$NF_est[,2:3],
  NM_est = full_data$NM_est[,2:3],
  # NF_ct = full_data$NF_ct,
  # NM_ct = full_data$NM_ct,
  r_dt = full_data$r_dt[,2:3],
  nc1e = full_data$nc1e,
  nf1e = full_data$nf1e,
  nm1e = full_data$nm1e,
  sc_cjs = ipm_data$s_cjs[62:92,4:5],
  sf_cjs = ipm_data$s_cjs[1:30, 4:5],
  sm_cjs = ipm_data$s_cjs[31:61,4:5],
  af_count = ipm_data$n_f_p_count[,4],
  am_count = ipm_data$n_m_p_count[,4]
  # l = full_data$l,
  # male = full_data$male,
  # female = full_data$female,
  # calf = full_data$calf,
  # herd = full_data$herd,
  # f = full_data$f,
  # nind = full_data$nind,
)

#Constants======================================================================
cnst <- list(
  NFmin = full_data$NFmin,
  NCmin = full_data$NCmin,
  NMmin = full_data$NMmin,
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
  nsc = nrow(dtf$sc_cjs),
  nsf = nrow(dtf$sf_cjs),
  nsm = nrow(dtf$sm_cjs),
  nn_fc = nrow(ipm_data$n_f_p_count),
  nn_mc = nrow(ipm_data$n_m_p_count),
  elk = full_data$elk_density,
  sc_cjs_t = ipm_data$s_cjs[62:92,1],
  sf_cjs_t = ipm_data$s_cjs[1:30, 1],
  sm_cjs_t = ipm_data$s_cjs[31:61,1],
  af_count_t = ipm_data$n_f_p_count[,1],
  am_count_t = ipm_data$n_m_p_count[,1],
  NC_est_t = full_data$NC_est[,1],
  NF_est_t = full_data$NF_est[,1],
  NM_est_t = full_data$NM_est[,1],
  r_dt_t = full_data$r_dt[,1]
)

#Initial values=================================================================

inits <- readRDS("data//the_ipm_inits.rds")

#Model setup====================================================================
mons <- c(
  "R", "SC", "SF", "SM", "NC", "NF", "NM", "R_Bvy", "R_Bvm", "R_B0_mean", 
  "R_Bpu", "R_Bdd", "R_B0_sd", "SC_B0_mean", "SC_B0_sd", "SC_Bvy", "SC_Bvm", 
  "SC_Bpu", "SC_Bdd", "NFaug", "NMaug", "SM_B0_mean", "SF_B0_mean", "SM_B0_sd",
  "SF_B0_sd", "sd_afcount", "sd_amcount",
  "NCaug", "Ntot"
)
# "PF", "PM",

# nimbleMCMC(code = model_code, data = dtf, constants = cnst, inits = inits,
#            monitors = mons, niter = 1000, nburnin = 100, nchains = 1)

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

block_variables <- function(cnf, variables, block_type){
  cnf$removeSamplers(variables)
  cnf$addSampler(variables, block_type)
}
sample_log <- function(cnf, variable){
  cnf$removeSampler(variable)
  cnf$addSampler(variable, "RW", control = list(log = T))
}
reflect_sample <- function(cnf, variable){
  cnf$removeSampler(variable)
  cnf$addSampler(variable, "RW", control = list(reflect = T))
}
reflect_sample <- function(cnf, variable){
  cnf$removeSampler(variable)
  cnf$addSampler(variable, "RW", control = list(reflect = T))
}
# block_variables(cnf, paste0("R_B0[", 2:34, "]"))
# map(c(
#   paste0("NCaug[", 2:32, "]"), 
#   paste0("NFaug[", 2:33, "]"), 
#   paste0("NMaug[", 2:33, "]")),
#   reflect_sample, cnf = cnf
# )
block_variables(
  cnf = cnf, 
  variables = paste0("NCaug[", 2:32, "]"), 
  block_type = "AF_slice"
)
block_variables(
  cnf = cnf, 
  variables = paste0("NFaug[", 2:33, "]"), 
  block_type = "AF_slice"
)
block_variables(
  cnf = cnf, 
  variables = paste0("NMaug[", 2:33, "]"), 
  block_type = "AF_slice"
)
# sample_log(cnf, paste0("NCaug[", 2:32, "]"))
# sample_log(cnf, paste0("NFaug[", 2:33, "]"))
# sample_log(cnf, paste0("NMaug[", 2:33, "]"))

# block_variables(cnf, paste0("SC_B0[", 2:32, "]"))
# block_variables(cnf, paste0("SF_B0[", 2:33, "]"))
# block_variables(cnf, paste0("SM_B0[", 2:33, "]"))

MCipm <- buildMCMC(cnf)
MCMCipm <- compileNimble(MCipm)
# c_ipm_mcmc$ipm_mcmc$run(10000)
# x <- c_ipm_mcmc$ipm_mcmc$mvSamples %>% as.matrix()

rslt <- runMCMC(
  MCMCipm, 
  niter = 100000, 
  nburnin = 0, 
  nchains = 5, 
  inits = inits,
  thin = 10
)
# c_ipm_mcmc$ipm_mcmc$run(100000, reset = FALSE)
# as.matrix(c_ipm_mcmc$ipm_mcmc$mvSamples) %>% mcmcplot()
mcmcplot(rslt)
