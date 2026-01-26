source("functions//calculate_lambda_functions.R")
tag <- "bst_lgst"
cg_cov_name <- case_when(
  tag == "bst_lgst" ~ "cg_logistic_growth",
  tag == "bst_mean" ~ "cg_mean_estimate",
  tag == "bst_odfw" ~ "cg_odfw_estimate",
  tag == "bst_rcns" ~ "cg_reconstruction"
)
load_file_name <- paste0("results//ipm_rslt_summ_lambda_prep_", tag, ".rds")
save_file_name <- paste0("results//full_lambda_tibble_", tag, ".rds")
tmls_file_name <- paste0("results//transition_matrix_list_", tag, ".rds")
covs_file_name <- paste0("data//tm_list_covariates_", tag, ".rds")
mtdt_file_name <- paste0("data//tm_list_metadata_", tag, ".rds")

rs <- readRDS(load_file_name) %>% 
  pivot_wider(
    id_cols = c(mcmc_step, chain, yr), 
    values_from = value, 
    names_from = param)
ipm_dt <- readRDS("data//ipm_data_30dec2025.rds")
cv <- tibble(
  cg = ipm_dt[[cg_cov_name]],
  dd = ipm_dt$nelk,
  wm = lag(ipm_dt$spei12),
  wt = ipm_dt$spei12
) %>%
  mutate(yr = 1:nrow(.))
saveRDS(cv, covs_file_name)

stp_vec <- rs %>% pull(mcmc_step) %>% unique()
yr_vec <- 2:36
pm_vec <- expand.grid(stpn = stp_vec, yrx = yr_vec)
if(tag == "bst_norm"){
  year_data <- rs %>% filter(!is.na(yr)) %>%
    mutate(stpn = mcmc_step, yrx = as.numeric(yr)) %>%
    select(stpn, yrx, ppb0o, ppb0y, ppb0p, ppyr, sn, snyr, sc, scyr, sf, sfyr, nf,
           nof, npf, nyf, nc, h1, h2, cges)
}else{
  year_data <- rs %>% filter(!is.na(yr)) %>%
    mutate(stpn = mcmc_step, yrx = as.numeric(yr)) %>%
    select(stpn, yrx, ppb0o, ppb0y, ppb0p, ppyr, sn, snyr, sc, scyr, sf, sfyr, nf,
           nof, npf, nyf, nc, h1, h2)
}
beta_data <- rs %>% filter(is.na(yr)) %>%
  mutate(stpn = mcmc_step) %>%
  select(
    stpn, 
                ppdd, ppwm,       ppb0y, ppb0p, ppb0o,
    snb0, sncg, sndd,       
    scb0, sccg,       scwm, scwt,
    sfb0,             sfwm, sfwt)
full_mtdt <- pm_vec %>%
  as_tibble() %>%
  left_join(year_data) %>%
  left_join(beta_data)

saveRDS(full_mtdt, mtdt_file_name)
plan(multisession, workers = 16)
tm_list <- future_pmap(.l = pm_vec, build_tm, rsdt = rs, cvdt = cv)
abundance_data <- full_mtdt %>%
  select(stpn, yrx, nc, nf, nof, npf, nyf)
Ndt <- pm_vec %>% as_tibble() %>%
  left_join(abundance_data) %>%
  # filter(stpn %in% 1:2) %>%
  select(nc, nf, nof, npf, nyf) %>% split(1:nrow(.))
tmna_list <- map2(.x = tm_list, .y = Ndt, .f = add_abundance_to_list)
saveRDS(tmna_list, file = tmls_file_name)
