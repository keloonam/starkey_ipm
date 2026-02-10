require(tidyverse); require(rjags)
tag <- "best_rcns"
folder <- "results//"
load_file <- paste0(folder, "fbipm_rslt_06feb2026_", tag, ".rds")
save_file <- paste0(folder, "ipm_rslt_summ_lambda_prep_", tag, ".rds")
rslt <- readRDS(load_file) %>%
  map(as.matrix) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(chain = c(
    rep("a", nrow(.)/3),
    rep("b", nrow(.)/3),
    rep("c", nrow(.)/3)
    )) %>%
  mutate(mcmc_step = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-2), names_to = "ipm_name") %>%
  mutate(param = case_when(
    grepl("SF_Byr", .$ipm_name) ~ "sfyr",
    grepl("SFB_cg", .$ipm_name) ~ "sfcg",
    grepl("SFB_dd", .$ipm_name) ~ "sfdd",
    grepl("SFB_wm", .$ipm_name) ~ "sfwm",
    grepl("SFB_wt", .$ipm_name) ~ "sfwt",
    grepl("SF_B0", .$ipm_name) ~ "sfb0",
    grepl("SF\\[", .$ipm_name) ~ "sf",
    grepl("SC_Byr", .$ipm_name) ~ "scyr",
    grepl("SCB_cg", .$ipm_name) ~ "sccg",
    grepl("SCB_dd", .$ipm_name) ~ "scdd",
    grepl("SCB_wm", .$ipm_name) ~ "scwm",
    grepl("SCB_wt", .$ipm_name) ~ "scwt",
    grepl("SC_B0", .$ipm_name) ~ "scb0",
    grepl("SC\\[", .$ipm_name) ~ "sc",
    grepl("SN_Byr", .$ipm_name) ~ "snyr",
    grepl("SNB_cg", .$ipm_name) ~ "sncg",
    grepl("SNB_dd", .$ipm_name) ~ "sndd",
    grepl("SNB_wt", .$ipm_name) ~ "snwt",
    grepl("SNB_wm", .$ipm_name) ~ "snwm",
    grepl("SN_B0", .$ipm_name) ~ "snb0",
    grepl("SN\\[", .$ipm_name) ~ "sn",
    grepl("P_Byr", .$ipm_name) ~ "ppyr",
    grepl("PB_cg", .$ipm_name) ~ "ppcg",
    grepl("PB_dd", .$ipm_name) ~ "ppdd",
    grepl("PB_wm", .$ipm_name) ~ "ppwm",
    grepl("P_B0\\[1", .$ipm_name) ~ "ppb0y",
    grepl("P_B0\\[2", .$ipm_name) ~ "ppb0p",
    grepl("P_B0\\[3", .$ipm_name) ~ "ppb0o",
    grepl("H\\[1,1", .$ipm_name) ~ "h1",
    grepl("H\\[2,1", .$ipm_name) ~ "h2",
    grepl("N_c\\[", .$ipm_name) ~ "nc",
    grepl("N_fo", .$ipm_name) ~ "nof",
    grepl("N_fp", .$ipm_name) ~ "npf",
    grepl("N_fy", .$ipm_name) ~ "nyf",
    grepl("N_f\\[", .$ipm_name) ~ "nf",
    grepl("cdrng\\[", .$ipm_name) ~ "cges",
    T ~ "remove"
  )) %>%
  filter(param != "remove")
params_without_years <- c(
  "ppb0", "snb0", "scb0", "sfb0", "ppwm", "ppdd", "ppwt", "ppcg", "ppb0",
  "snwt", "snwm", "sndd", "sncg", "scdd", "sccg", "scwm", "scwt", "sfcg", 
  "sfdd", "sfwm", "sfwt")
pull_yr <- function(dtx, x){
  dtx %>%
    separate(ipm_name, into = c(NA, "yr"), x) %>%
    separate(yr, into = c("yr", NA), -1) %>%
    pull(yr) %>%
    return()
}
rslt2 <- rslt %>% 
  mutate(has_year = case_when(
    param %in% c(params_without_years) ~ F,
    T ~ T
  )) %>%
  mutate(yr = case_when(
    has_year == F ~ NA,
    param %in% c("scyr", "sfyr", "snyr") ~ pull_yr(., 7),
    param %in% c("h1", "h2", "ppyr", "cges") ~ pull_yr(., 6),
    param %in% c("nof", "nyf", "npf") ~ pull_yr(., 5),
    param %in% c("nf", "nc") ~ pull_yr(., 4),
    param %in% c("sc", "sf", "sn") ~ pull_yr(., 3),
    param %in% c("pp") ~ pull_yr(., 2)
  )) %>%
  select(mcmc_step, chain, param, yr, value)
saveRDS(rslt2, save_file)
rm(list = ls())
