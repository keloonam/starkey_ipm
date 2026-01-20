#Environment====================================================================
require(tidyverse); source("functions//lambda_analysis_functions.R")

#Shared Objects=================================================================
all_params <- c(
  "nfp", "ncp", "sh1", "sh2", "cg", "dd", "wt", "wm", "cm", "dm",
  "ppb0", "ppcg", "ppwm", "ppyr", "ppdd",
  "snb0", "sncg", "snwt", "snyr", "sndd",
  "scb0", "sccg", "scwt", "scwm", "scdd", "scyr",
  "sfb0", "sfcg", "sfwt", "sfwm", "sfdd", "sfyr"
  )
params_oi <- c("nfp", "ncp", "sh1", "sh2", "cg", "dd", "wt", "wm", "cm", "dm",
               "ppyr", "snyr", "scyr", "sfyr")
mdt <- readRDS("data//tm_list_covariates.rds") %>%
  mutate(cm = lag(cg), dm = lag(dd)) %>%
  mutate(yrx = yr) %>% select(-yr)
fdt <- readRDS("results//full_lambda_tibble.rds") %>% 
  as_tibble() %>% 
  left_join(mdt) %>%
  mutate(sh1 = 1-h1, sh2 = 1-h2) %>%
  mutate(nfp = nf/(nc+nf), ncp = nc/(nc+nf))

Leq <- expression(
  # Numerator
  # Reproduction
  (0.5 *
  (1 /(1 + exp(-(ppb0 + ppyr + ppcg*cm + ppdd*dm + ppwm*wm)))) * 
  (1 /(1 + exp(-(snb0 + snyr + sncg*cg + sndd*dd + snwt*wt)))) * (
  nfp * sh2 /(1 + exp(-(sfb0 + sfyr + sfcg*cg + sfdd*dd + sfwm*wm + sfwt*wt))) +
  ncp * sh1 /(1 + exp(-(scb0 + scyr + sccg*cg + scdd*dd + scwm*wm + scwt*wt))))+
  # Survival
  nfp * sh2 /(1 + exp(-(sfb0 + sfyr + sfcg*cg + sfdd*dd + sfwm*wm + sfwt*wt))) +
  ncp * sh1 /(1 + exp(-(scb0 + scyr + sccg*cg + scdd*dd + scwm*wm + scwt*wt))))/
  # Denominator
  (nfp + ncp)
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
sensitivities %>% saveRDS("results//lambda_sensitivity.rds")

list(wvcv = elasticities_wvcv, solo = elasticities_solo) %>%
  saveRDS("results//lambda_elasticity.rds")

#Annual Variability=============================================================
plan(multisession, workers = 12)
dmdt <- future_map(.x = unique(fdt$stpn), .f = get_dfs_mns, fdt = fdt) %>%
  bind_rows()
dfs <- array(
  data = NA, 
  dim = c(length(unique(fdt$stpn)), length(unique(dmdt$yrx)), length(params_oi))
)
dimnames(dfs)[[3]] <- params_oi
dfs[,,"nfp"] <- pull_diff_mat(dmdt, "dnf")
dfs[,,"ncp"] <- pull_diff_mat(dmdt, "dnc")
dfs[,,"sh1"] <- pull_diff_mat(dmdt, "dh1")
dfs[,,"sh2"] <- pull_diff_mat(dmdt, "dh2")
dfs[,,"cg"] <- pull_diff_mat(dmdt, "dcg")
dfs[,,"dd"] <- pull_diff_mat(dmdt, "ddd")
dfs[,,"wt"] <- pull_diff_mat(dmdt, "dwt")
dfs[,,"wm"] <- pull_diff_mat(dmdt, "dwm")
dfs[,,"cm"] <- pull_diff_mat(dmdt, "dcm")
dfs[,,"dm"] <- pull_diff_mat(dmdt, "ddm")
dfs[,,"ppyr"] <- pull_diff_mat(dmdt, "dpp")
dfs[,,"snyr"] <- pull_diff_mat(dmdt, "dsn")
dfs[,,"scyr"] <- pull_diff_mat(dmdt, "dsc")
dfs[,,"sfyr"] <- pull_diff_mat(dmdt, "dsf")

mns <- map(.x = all_params, .f = pull_mn_mat, dtx = dmdt)
names(mns) <- all_params

sns <- array(NA, dim = dim(dfs))
dimnames(sns)[[3]] <- params_oi
sns[,,"nfp"] <- calc_deriv_ann(param = "nfp", eq = Leq, mu = mns)
sns[,,"ncp"] <- calc_deriv_ann(param = "ncp", eq = Leq, mu = mns)
sns[,,"sh1"] <- calc_deriv_ann(param = "sh1", eq = Leq, mu = mns)
sns[,,"sh2"] <- calc_deriv_ann(param = "sh2", eq = Leq, mu = mns)
sns[,,"cg"] <- calc_deriv_ann(param = "cg", eq = Leq, mu = mns)
sns[,,"dd"] <- calc_deriv_ann(param = "dd", eq = Leq, mu = mns)
sns[,,"wt"] <- calc_deriv_ann(param = "wt", eq = Leq, mu = mns)
sns[,,"wm"] <- calc_deriv_ann(param = "wm", eq = Leq, mu = mns)
sns[,,"cm"] <- calc_deriv_ann(param = "cm", eq = Leq, mu = mns)
sns[,,"dm"] <- calc_deriv_ann(param = "dm", eq = Leq, mu = mns)
sns[,,"ppyr"] <- calc_deriv_ann(param = "ppyr", eq = Leq, mu = mns)
sns[,,"snyr"] <- calc_deriv_ann(param = "snyr", eq = Leq, mu = mns)
sns[,,"scyr"] <- calc_deriv_ann(param = "scyr", eq = Leq, mu = mns)
sns[,,"sfyr"] <- calc_deriv_ann(param = "sfyr", eq = Leq, mu = mns)

annual_effect <- sns * dfs
saveRDS(annual_effect, file = "results//annual_L_effects.rds")
