source("functions//calculate_lambda_functions.R")

rs <- readRDS("results//ipm_rslt_summ_lambda_prep.rds") %>% 
  pivot_wider(
    id_cols = c(mcmc_step, yr), 
    values_from = value, 
    names_from = param)
ipm_dt <- readRDS("data//ipm_data_17dec2025.rds")
cv <- tibble(
  cg = ipm_dt$cdens,
  dd = ipm_dt$nelk,
  wm = lag(ipm_dt$spei12),
  wt = ipm_dt$spei12
) %>%
  mutate(yr = 1:nrow(.))
saveRDS(cv, "data//tm_list_covariates.rds")

stp_vec <- rs %>% pull(mcmc_step) %>% unique()
yr_vec <- 2:36
pm_vec <- expand.grid(stpn = stp_vec, yrx = yr_vec)
year_data <- rs %>% filter(!is.na(yr)) %>%
  mutate(stpn = mcmc_step, yrx = as.numeric(yr)) %>%
  select(stpn, yrx, pp, ppyr, sn, snyr, sc, scyr, sf, sfyr, nf, nc, h1, h2)
beta_data <- rs %>% filter(is.na(yr)) %>%
  mutate(stpn = mcmc_step) %>%
  select(
    stpn, 
    ppb0, ppcg, ppdd, ppwm, 
    snb0, sncg, sndd,       snwt,
    scb0, sccg, scdd, scwm, scwt,
    sfb0, sfcg, sfdd, sfwm, sfwt)
full_mtdt <- pm_vec %>%
  as_tibble() %>%
  left_join(year_data) %>%
  left_join(beta_data)
saveRDS(full_mtdt, "data//tm_list_metadata.rds")
tm_list <- pmap(.l = pm_vec, build_tm, rsdt = rs, cvdt = cv)
tm_list <- readRDS("results//transition_matrix_list.rds")
abundance_data <- full_mtdt %>%
  select(stpn, yrx, nc, nf)
Ndt <- pm_vec %>% as_tibble() %>%
  left_join(abundance_data) %>%
  # filter(stpn %in% 1:2) %>%
  select(nc, nf) %>% split(1:nrow(.))
tmna_list <- map2(.x = tm_list, .y = Ndt, .f = add_abundance_to_list)
saveRDS(tmna_list, file = "results//transition_matrix_list.rds")

# tmna_list <- readRDS("results//transition_matrix_list.rds")
plan(multisession, workers = 12)
lambda_results <- future_map(
  .x = tmna_list,
  .f = calculate_lambdas
) %>%
  bind_rows()
saveRDS("results//lambda_rslt.rds")
# lambda_results <- readRDS("results//lambda_rslt.rds")
cbind(lambda_results, full_mtdt) %>%
  saveRDS("results//full_lambda_tibble.rds")
