#Environment====================================================================
require(tidyverse); source("functions//lambda_analysis_functions.R")
tag <- "best_rcns"
lam_tib_name <- paste0("results//full_lambda_tibble_", tag, ".rds")
met_dat_name <- paste0("data//tm_list_covariates_", tag, ".rds")
sens_save_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
elas_save_name <- paste0("results//lambda_elasticity_", tag, ".rds")

use_supported_only <- TRUE
#Shared Objects=================================================================
all_params <- c(
  "pr2", "prc", "pry", "prp", "pro", 
  "sh1", "sh2", 
  "cg", "dd", "wt", "wm",
  "ppwm", "ppyr", "ppdd", "ppb0o", "ppb0y", "ppb0p", 
  "snb0", "snwm", "sncg", "sndd", "snyr",
  "scb0", "scwt", "scwm", "scyr", "sccg", 
  "sfb0", "sfwm", "sfyr", "sfdd"
  )
params_oi <- c("pr2", "prc", "pry", "prp", "pro", 
               "sh1", "sh2", 
               "cg", "dd", "wt", "wm", 
               "ppyr", "snyr", "scyr", "sfyr")
mdt <- readRDS(met_dat_name) %>%
  mutate(yrx = yr) %>% select(-yr)
fdt <- readRDS(lam_tib_name) %>% 
  as_tibble() %>% 
  left_join(mdt) %>%
  mutate(sh1 = 1-h1, sh2 = 1-h2) %>%
  mutate(
    pry = nyf / (nf + nc),
    prp = npf / (nf + nc),
    pro = nof / (nf + nc),
    pr2 = (nf - nyf - npf - nof)  / (nf + nc),
    prc = nc  / (nf + nc)
  )
if("cges" %in% names(fdt)){
  fdt <- fdt %>%
    mutate(cg = cges)
}

Leq <- expression(
  # Numerator
  # Reproduction
  (0.5 *
  ((1 /(1 + exp(-(ppb0y + ppyr + ppdd*dd + ppwm*wm)))) * pry +
   (1 /(1 + exp(-(ppb0p + ppyr + ppdd*dd + ppwm*wm)))) * prp +
   (1 /(1 + exp(-(ppb0o + ppyr + ppdd*dd + ppwm*wm)))) * pro) *
  (1 /(1 + exp(-(snb0 + snyr + sncg*cg + sndd*dd + snwm*wm)))) * 
  ((pr2+pry+prp+pro) * sh2 / 
     (1 + exp(-(sfb0 + sfyr + sfwm*wm + sfdd*dd))) +
  prc * sh1 / 
    (1 + exp(-(scb0 + scyr + sccg*cg + scwm*wm + scwt*wt)))) +
  # Survival
  (pr2+pry+prp+pro) * sh2 / 
    (1 + exp(-(sfb0 + sfyr + sfwm*wm + sfdd*dd))) +
  prc * sh1 /(1 + exp(-(scb0 + scyr + sccg*cg + scwm*wm + scwt*wt))))/
  # Denominator
  ((pr2+pry+prp+pro) + prc)
)

#Total variability==============================================================
mu <- map(.x = all_params, .f = pull_param_vec, fdt)
names(mu) <- all_params

sensitivities <- map(
  .x = params_oi, 
  .f = calculate_derivative, 
  eq = Leq, 
  mu = mu) %>%
  bind_rows()
require(furrr)
plan(multisession, workers = 12)
elasticities_wvcv <- future_map(
  .x = unique(fdt$stpn),
  .f = calculate_elasticity_with_covariances,
  fdt = fdt,
  sensitivities = sensitivities
) %>%
  bind_rows() 
plan(multisession, workers = 12)
elasticities_solo <- future_map(
  .x = unique(fdt$stpn),
  .f = calculate_elasticity_without_covariances,
  fdt = fdt,
  sensitivities = sensitivities
) %>%
  bind_rows() 
sensitivities %>% saveRDS(sens_save_name)

list(wvcv = elasticities_wvcv, solo = elasticities_solo) %>%
  saveRDS(elas_save_name)
rm(list = ls())
#Annual Variability=============================================================
# plan(multisession, workers = 12)
# dmdt <- future_map(.x = unique(fdt$stpn), .f = get_dfs_mns, fdt = fdt) %>%
#   bind_rows()
# dfs <- array(
#   data = NA, 
#   dim = c(length(unique(fdt$stpn)), length(unique(dmdt$yrx)), length(params_oi))
# )
# dimnames(dfs)[[3]] <- params_oi
# dfs[,,"nfp"] <- pull_diff_mat(dmdt, "dnf")
# dfs[,,"ncp"] <- pull_diff_mat(dmdt, "dnc")
# dfs[,,"sh1"] <- pull_diff_mat(dmdt, "dh1")
# dfs[,,"sh2"] <- pull_diff_mat(dmdt, "dh2")
# dfs[,,"cg"] <- pull_diff_mat(dmdt, "dcg")
# dfs[,,"dd"] <- pull_diff_mat(dmdt, "ddd")
# dfs[,,"wt"] <- pull_diff_mat(dmdt, "dwt")
# dfs[,,"wm"] <- pull_diff_mat(dmdt, "dwm")
# dfs[,,"cm"] <- pull_diff_mat(dmdt, "dcm")
# dfs[,,"dm"] <- pull_diff_mat(dmdt, "ddm")
# dfs[,,"ppyr"] <- pull_diff_mat(dmdt, "dpp")
# dfs[,,"snyr"] <- pull_diff_mat(dmdt, "dsn")
# dfs[,,"scyr"] <- pull_diff_mat(dmdt, "dsc")
# dfs[,,"sfyr"] <- pull_diff_mat(dmdt, "dsf")
# 
# mns <- map(.x = all_params, .f = pull_mn_mat, dtx = dmdt)
# names(mns) <- all_params
# 
# sns <- array(NA, dim = dim(dfs))
# dimnames(sns)[[3]] <- params_oi
# sns[,,"nfp"] <- calc_deriv_ann(param = "nfp", eq = Leq, mu = mns)
# sns[,,"ncp"] <- calc_deriv_ann(param = "ncp", eq = Leq, mu = mns)
# sns[,,"sh1"] <- calc_deriv_ann(param = "sh1", eq = Leq, mu = mns)
# sns[,,"sh2"] <- calc_deriv_ann(param = "sh2", eq = Leq, mu = mns)
# sns[,,"cg"] <- calc_deriv_ann(param = "cg", eq = Leq, mu = mns)
# sns[,,"dd"] <- calc_deriv_ann(param = "dd", eq = Leq, mu = mns)
# sns[,,"wt"] <- calc_deriv_ann(param = "wt", eq = Leq, mu = mns)
# sns[,,"wm"] <- calc_deriv_ann(param = "wm", eq = Leq, mu = mns)
# sns[,,"cm"] <- calc_deriv_ann(param = "cm", eq = Leq, mu = mns)
# sns[,,"dm"] <- calc_deriv_ann(param = "dm", eq = Leq, mu = mns)
# sns[,,"ppyr"] <- calc_deriv_ann(param = "ppyr", eq = Leq, mu = mns)
# sns[,,"snyr"] <- calc_deriv_ann(param = "snyr", eq = Leq, mu = mns)
# sns[,,"scyr"] <- calc_deriv_ann(param = "scyr", eq = Leq, mu = mns)
# sns[,,"sfyr"] <- calc_deriv_ann(param = "sfyr", eq = Leq, mu = mns)
# 
# annual_effect <- sns * dfs
# saveRDS(annual_effect, file = "results//annual_L_effects.rds")
