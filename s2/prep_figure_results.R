#Header=========================================================================
# Kenneth Loonam
# July 2025

#Environment====================================================================
# Packages
require(dplyr); require(rjags); require(readr); require(tidyr); require(purrr)
require(lubridate); require(ggplot2)
# Functions
source("s2//plotting_functions.R")

# Results to load
ipm_rslt <- readRDS("s2//results//ipmrs_27sep2025_spei12.rds") %>%
  map(as_tibble) %>%
  bind_rows()
ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")
preg <- readRDS("results//pregnancy_analysis_13oct2024.rds") %>%
  map(as_tibble) %>%
  bind_rows()

#Demographic data===============================================================
n_ca <- pull_demo(ipm_rslt, pat = "N_c", yrs = 1988:2023)
n_af <- pull_demo(ipm_rslt, pat = "N_f", yrs = 1988:2023)
n_am <- pull_demo(ipm_rslt, pat = "N_m", yrs = 1988:2023)
n_tot <- pull_demo(ipm_rslt, pat = "N_tot", yrs = 1988:2023)
r_dat <- pull_demo(ipm_rslt, pat = "R\\[", yrs = 1989:2023)
sf_dat <- pull_demo(ipm_rslt, pat = "survival_af", yrs = 1989:2023)
sm_dat <- pull_demo(ipm_rslt, pat = "survival_am", yrs = 1989:2023)
sc_dat <- pull_demo(ipm_rslt, pat = "survival_ca", yrs = 1989:2023)

#Covariate data=================================================================
cov_dat <- tibble(
  covariate = c(
    rep("spei_t", 36), 
    rep("cg_dns", 36), 
    rep("af_dns", 36)),
  value = c(
    ipm_data$spei12, 
    ipm_data$cdens, 
    ipm_data$nelk),
  year = c(1988:2023, 1988:2023, 1988:2023)
) %>%
  pivot_wider(names_from = covariate, values_from = value) %>%
  mutate(spei_m = lag(spei_t, 1))

cougar_density <-  c(
  0.14336299, 0.19115066, 0.23893832, 0.14336299, 0.23893832, 0.38230132,
  0.40671006, 0.67690686, 0.92058393, 1.21493542, 1.38220093, 1.60460137, 
  1.81466363, 1.99925826, 2.01422234, 2.01290939, 2.15601154, 2.10612837, 
  2.12870926, 1.98718091, 1.99400648, 1.99400648, 1.99689322, 1.94568843, 
  2.13999536, 1.89817030, 1.87585895, 1.97222552, 1.89214465, 1.81652433, 
  1.89162295, 1.90606533, 2.03551633, 2.00821406, 1.90843904, 1.62328685)
cougar_dat <- tibble(
  year = 1988:2023,
  cougar_density = cougar_density
)
starkey_area_km2 <- 77.382817
n_af_cov <- n_af %>% mutate(year = yr) %>% 
  mutate(mean_density = mci) %>%
  select(mean_density, year)
spei_real_scale <- read_csv("data//climate//spei.csv") %>%
  mutate(dt = my(DATA)) %>%
  mutate(yr = year(dt)) %>%
  mutate(mn = month(dt)) %>%
  mutate(spei_1m = SPEI_1) %>%
  mutate(spei_3m = SPEI_3) %>%
  mutate(spei_6m = SPEI_6) %>%
  mutate(spei_12m = SPEI_12) %>%
  mutate(spei_24m = SPEI_24) %>%
  mutate(spei_48m = SPEI_48) %>%
  filter(yr > 1987 & yr < 2024) %>%
  filter(mn == 9) %>%
  select(yr, spei_1m, spei_3m, spei_6m, spei_12m, spei_24m, spei_48m) %>%
  arrange(yr) %>%
  pull(spei_12m)
spei_dat <- tibble(
  year = 1988:2023,
  spei = spei_real_scale
)

cov_dat_real <- full_join(cougar_dat, spei_dat) %>% full_join(n_af_cov) %>%
  mutate(
    cg = cougar_density,
    wm = lag(spei),
    wt = spei,
    ed = mean_density) %>%
  select(year, cg, wm, wt, ed)

#Partial residual point data====================================================
# Recruitment
rbdt <- ipm_rslt %>%
  mutate(id = 1:nrow(.)) %>%
  select(c("R_cg", "R_dd", "R_wm", "R_wt", "id"))
r_point_dat <- ipm_rslt %>% select(grep("R\\[", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "r") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(rbdt) %>% left_join(cov_dat) %>%
  mutate(b0r = logit(r)) %>%
  mutate(cg_resid = expit(b0r - R_dd*af_dns - R_wm*spei_m - R_wt*spei_t)) %>%
  mutate(wt_resid = expit(b0r - R_dd*af_dns - R_wm*spei_m - R_cg*cg_dns)) %>%
  mutate(wm_resid = expit(b0r - R_dd*af_dns - R_cg*cg_dns - R_wt*spei_t)) %>%
  mutate(ed_resid = expit(b0r - R_cg*cg_dns - R_wm*spei_m - R_wt*spei_t)) %>%
  select(year, cg_resid, wt_resid, wm_resid, ed_resid, r) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "param", values_to = "val") %>%
  group_by(year, param) %>%
  summarise(
    lci = quantile(val, 0.025),
    mci = quantile(val, 0.500),
    uci = quantile(val, 0.975)
  ) %>%
  left_join(cov_dat_real)

# Survival
sbdt <- ipm_rslt %>%
  mutate(id = 1:nrow(.)) %>%
  select(c("S_cg", "S_dd", "S_wm", "S_wt", "id"))
s_point_dat <- ipm_rslt %>% select(grep("survival_ca", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "s") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(sbdt) %>% left_join(cov_dat) %>%
  mutate(b0s = logit(s)) %>%
  mutate(cg_resid = expit(b0s - S_dd*af_dns - S_wm*spei_m - S_wt*spei_t)) %>%
  mutate(wt_resid = expit(b0s - S_dd*af_dns - S_wm*spei_m - S_cg*cg_dns)) %>%
  mutate(wm_resid = expit(b0s - S_dd*af_dns - S_cg*cg_dns - S_wt*spei_t)) %>%
  mutate(ed_resid = expit(b0s - S_cg*cg_dns - S_wm*spei_m - S_wt*spei_t)) %>%
  select(year, cg_resid, wt_resid, wm_resid, ed_resid, s) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "param", values_to = "val") %>%
  group_by(year, param) %>%
  summarise(
    lci = quantile(val, 0.025),
    mci = quantile(val, 0.500),
    uci = quantile(val, 0.975)
  ) %>%
  left_join(cov_dat_real)

# Pregnancy
pbdt <- preg %>%
  mutate(id = 1:nrow(.)) %>%
  select(c("bden_lac_prm", "bpdi_lac_prm", "id"))
p_point_dat <- preg %>% select(grep("a_preg_lac_prime", names(.))) %>%
  set_names(as.character(1989:2021)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "p") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(pbdt) %>% left_join(cov_dat) %>%
  mutate(b0p = logit(p)) %>%
  mutate(wt_resid = expit(b0p - bden_lac_prm*af_dns)) %>%
  mutate(ed_resid = expit(b0p - bpdi_lac_prm*spei_t)) %>%
  select(year, wt_resid, ed_resid, p) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "param", values_to = "val") %>%
  group_by(year, param) %>%
  summarise(
    lci = quantile(val, 0.025),
    mci = quantile(val, 0.500),
    uci = quantile(val, 0.975)
  ) %>%
  left_join(cov_dat_real)

# Lambda
sca_tmp <- ipm_rslt %>% select(grep("survival_ca", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "s") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(sbdt) %>% left_join(cov_dat) %>%
  mutate(b0s = logit(s)) %>%
  mutate(cg_resid_s = expit(b0s - S_dd*af_dns - S_wm*spei_m - S_wt*spei_t)) %>%
  mutate(wt_resid_s = expit(b0s - S_dd*af_dns - S_wm*spei_m - S_cg*cg_dns)) %>%
  mutate(wm_resid_s = expit(b0s - S_dd*af_dns - S_cg*cg_dns - S_wt*spei_t)) %>%
  mutate(ed_resid_s = expit(b0s - S_cg*cg_dns - S_wm*spei_m - S_wt*spei_t))
rdt_tmp <- ipm_rslt %>% select(grep("R\\[", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "r") %>%
  mutate(year = as.numeric(year)) %>%
  left_join(rbdt) %>% left_join(cov_dat) %>%
  mutate(b0r = logit(r)) %>%
  mutate(cg_resid_r = expit(b0r - R_dd*af_dns - R_wm*spei_m - R_wt*spei_t)) %>%
  mutate(wt_resid_r = expit(b0r - R_dd*af_dns - R_wm*spei_m - R_cg*cg_dns)) %>%
  mutate(wm_resid_r = expit(b0r - R_dd*af_dns - R_cg*cg_dns - R_wt*spei_t)) %>%
  mutate(ed_resid_r = expit(b0r - R_cg*cg_dns - R_wm*spei_m - R_wt*spei_t))
afs_tmp <- ipm_rslt %>% select(grep("survival_af", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "afs") %>%
  mutate(year = as.numeric(year))
afn_tmp <- ipm_rslt %>% select(grep("N_f", names(.))) %>%
  set_names(as.character(1988:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "afn") %>%
  mutate(year = as.numeric(year))
can_tmp <- ipm_rslt %>% select(grep("N_c", names(.))) %>%
  set_names(as.character(1988:2023)) %>%
  mutate(id = 1:nrow(.)) %>%
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "year", values_to = "can") %>%
  mutate(year = as.numeric(year))

ddt_tmp <- full_join(sca_tmp, rdt_tmp) %>% full_join(afs_tmp) %>%
  left_join(afn_tmp) %>% left_join(can_tmp) %>%
  filter(id %in% sample(unique(.$id), 5000)) %>%
  rowwise() %>%
  mutate(fl_nu = calc_e_lam(fs = afs, r = r, cs = s)) %>%
  mutate(tl_nu = calc_t_lam(fs = afs, r = r, cs = s, nc = can, nf = afn)) %>%
  mutate(fl_cg = calc_e_lam(fs = afs, r = cg_resid_r, cs = cg_resid_s)) %>%
  mutate(fl_wm = calc_e_lam(fs = afs, r = wm_resid_r, cs = wm_resid_s)) %>%
  mutate(fl_wt = calc_e_lam(fs = afs, r = wt_resid_r, cs = wt_resid_s)) %>%
  mutate(fl_ed = calc_e_lam(fs = afs, r = ed_resid_r, cs = wt_resid_s)) %>%
  mutate(tl_cg = calc_t_lam(fs = afs, r = cg_resid_r, cs = cg_resid_s, 
                            nc = can, nf = afn)) %>%
  mutate(tl_wm = calc_t_lam(fs = afs, r = wm_resid_r, cs = wm_resid_s, 
                            nc = can, nf = afn)) %>%
  mutate(tl_wt = calc_t_lam(fs = afs, r = wt_resid_r, cs = wt_resid_s, 
                            nc = can, nf = afn)) %>%
  mutate(tl_ed = calc_t_lam(fs = afs, r = ed_resid_r, cs = wt_resid_s, 
                            nc = can, nf = afn)) %>%
  ungroup() %>%
  select(year, fl_nu, fl_cg, fl_wm, fl_wt, fl_ed, 
               tl_nu, tl_cg, tl_wm, tl_wt, tl_ed) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "param", values_to = "val") %>%
  group_by(year, param) %>%
  summarise(
    lci = quantile(val, 0.025), 
    mci = quantile(val, 0.5), 
    uci = quantile(val, 0.975)) %>%
  left_join(cov_dat_real)

list(
  lambda = ddt_tmp, 
  recruitment = r_point_dat, 
  survival = s_point_dat, 
  pregnancy = p_point_dat
  ) %>% saveRDS(file = "results//figure_partial_residuals.rds")

rm(list = ls())

#Line data======================================================================
#Load data======================================================================
# Functions
source("s2//plotting_functions.R")

# Results to load
ipm_rslt <- readRDS("s2//results//ipmrs_27sep2025_spei12.rds") %>%
  map(as_tibble) %>%
  bind_rows()
ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")
preg <- readRDS("results//pregnancy_analysis_13oct2024.rds") %>%
  map(as_tibble) %>%
  bind_rows()

# Beta results
cov_summaries <- ipm_rslt %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_dd", "R_wm", "R_wt", 
                          "S_cg", "S_dd", "S_wm", "S_wt")) %>%
  mutate(direction = value > 0) %>%
  group_by(parameter) %>%
  summarise(
    lcri = quantile(value, .025),
    medi = quantile(value, 0.50),
    ucri = quantile(value, .975),
    pofd = mean(direction)
  )



rb_post <- ipm_rslt %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_B0", "R_dd", "R_wm", "R_wt")) %>%
  mutate(id = sort(rep(1:(nrow(.)/5), 5))) %>%
  pivot_wider(id_cols = id, names_from = parameter)
sb_post <- ipm_rslt %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("S_cg", "S_C_B0_ps", "S_dd", "S_wm", "S_wt")) %>%
  mutate(id = sort(rep(1:(nrow(.)/5), 5))) %>%
  pivot_wider(id_cols = id, names_from = parameter) %>%
  mutate(S_B0 = logit(S_C_B0_ps)) %>%
  select(S_B0, S_cg, S_dd, S_wm, S_wt)
pb_post <- preg %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("b_lac_prm_mean", "b0_mean", 
                          "bden_lac_prm", "bpdi_lac_prm")) %>%
  mutate(id = sort(rep(1:(nrow(.)/4), 4))) %>%
  pivot_wider(id_cols = id, names_from = parameter)

# Marginal Plot Data
cov_vals <- seq(-3, 3, length.out = 1000)
### Recruitment ###
# Puma
r_cg_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = rb_post$R_B0, 
  b1 = rb_post$R_cg)
# SPEI t-1
r_wm_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = rb_post$R_B0, 
  b1 = rb_post$R_wm)
# SPEI
r_wt_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = rb_post$R_B0, 
  b1 = rb_post$R_wt)
# Density dependence
r_dd_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = rb_post$R_B0, 
  b1 = rb_post$R_dd)

### Survival ###
# Puma
s_cg_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = sb_post$S_B0, 
  b1 = sb_post$S_cg)
# SPEI t-1
s_wm_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = sb_post$S_B0, 
  b1 = sb_post$S_wm)
# SPEI
s_wt_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = sb_post$S_B0, 
  b1 = sb_post$S_wt)
# Density dependence
s_dd_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = sb_post$S_B0, 
  b1 = sb_post$S_dd)

### Lambda ###
# Puma
cg_lam_dt <- sim_lam(
  fdt = ipm_rslt, 
  cov = ipm_data$cdens, 
  length_cov = 500,
  r_cov_name = "R_cg",
  s_cov_name = "S_cg")
# SPEI t-1
wm_lam_dt <- sim_lam(
  fdt = ipm_rslt, 
  cov = ipm_data$spei12, 
  length_cov = 500,
  r_cov_name = "R_wm",
  s_cov_name = "S_wm")
# SPEI
wt_lam_dt <- sim_lam(
  fdt = ipm_rslt, 
  cov = ipm_data$spei12, 
  length_cov = 500,
  r_cov_name = "R_wt",
  s_cov_name = "S_wt")
# Density dependence
dd_lam_dt <- sim_lam(
  fdt = ipm_rslt, 
  cov = ipm_data$nelk, 
  length_cov = 500,
  r_cov_name = "R_dd",
  s_cov_name = "S_dd")
### Pregnancy ###
# SPEI
p_wt_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = pb_post$b0_mean + pb_post$b_lac_prm_mean, 
  b1 = pb_post$bpdi_lac_prm)
# Density dependence
p_dd_line <- calc_marg_line(
  cov = cov_vals, 
  b0 = pb_post$b0_mean + pb_post$b_lac_prm_mean, 
  b1 = pb_post$bden_lac_prm)


n_af <- pull_demo(ipm_rslt, pat = "N_f", yrs = 1988:2023)


cougar_density <-  c(
  0.14336299, 0.19115066, 0.23893832, 0.14336299, 0.23893832, 0.38230132,
  0.40671006, 0.67690686, 0.92058393, 1.21493542, 1.38220093, 1.60460137, 
  1.81466363, 1.99925826, 2.01422234, 2.01290939, 2.15601154, 2.10612837, 
  2.12870926, 1.98718091, 1.99400648, 1.99400648, 1.99689322, 1.94568843, 
  2.13999536, 1.89817030, 1.87585895, 1.97222552, 1.89214465, 1.81652433, 
  1.89162295, 1.90606533, 2.03551633, 2.00821406, 1.90843904, 1.62328685)
cougar_dat <- tibble(
  year = 1988:2023,
  cougar_density = cougar_density
)
starkey_area_km2 <- 77.382817
n_af_cov <- n_af %>% mutate(year = yr) %>% 
  mutate(mean_density = mci) %>%
  select(mean_density, year)
spei_real_scale <- read_csv("data//climate//spei.csv") %>%
  mutate(dt = my(DATA)) %>%
  mutate(yr = year(dt)) %>%
  mutate(mn = month(dt)) %>%
  mutate(spei_1m = SPEI_1) %>%
  mutate(spei_3m = SPEI_3) %>%
  mutate(spei_6m = SPEI_6) %>%
  mutate(spei_12m = SPEI_12) %>%
  mutate(spei_24m = SPEI_24) %>%
  mutate(spei_48m = SPEI_48) %>%
  filter(yr > 1987 & yr < 2024) %>%
  filter(mn == 9) %>%
  select(yr, spei_1m, spei_3m, spei_6m, spei_12m, spei_24m, spei_48m) %>%
  arrange(yr) %>%
  pull(spei_12m)
spei_dat <- tibble(
  year = 1988:2023,
  spei = spei_real_scale
)

cov_dat_real <- full_join(cougar_dat, spei_dat) %>% full_join(n_af_cov) %>%
  mutate(
    cg = cougar_density,
    wm = lag(spei),
    wt = spei,
    ed = mean_density) %>%
  select(year, cg, wm, wt, ed)



cgmn <- mean(cov_dat_real$cg)
cgsd <-   sd(cov_dat_real$cg)
spmn <- mean(cov_dat_real$wt)
spsd <-   sd(cov_dat_real$wt)
edmn <- mean(cov_dat_real$ed)
edsd <-   sd(cov_dat_real$ed) 


r_cg_ln <- r_cg_line %>%
  mutate(val = val * cgsd + cgmn)
r_pm_ln <- r_wm_line %>%
  mutate(val = val * spsd + spmn)
r_pt_ln <- r_wt_line %>%
  mutate(val = val * spsd + spmn)
r_ed_ln <- r_dd_line %>%
  mutate(val = val * edsd + edmn)
s_cg_ln <- s_cg_line %>%
  mutate(val = val * cgsd + cgmn)
s_pm_ln <- s_wm_line %>%
  mutate(val = val * spsd + spmn)
s_pt_ln <- s_wt_line %>%
  mutate(val = val * spsd + spmn)
s_ed_ln <- s_dd_line %>%
  mutate(val = val * edsd + edmn)
lam_cg_line <- cg_lam_dt$lam_df %>%
  mutate(val = val * cgsd + cgmn)
lam_pm_line <- wm_lam_dt$lam_df %>%
  mutate(val = val * spsd + spmn)
lam_pt_line <- wt_lam_dt$lam_df %>%
  mutate(val = val * spsd + spmn)
lam_ed_line <- dd_lam_dt$lam_df %>%
  mutate(val = val * edsd + edmn)
p_wt_ln <- p_wt_line %>%
  mutate(val = val * spsd + spmn)
p_ed_ln <- p_dd_line %>%
  mutate(val = val * edsd + edmn)

list(
  rcg_ln = r_cg_ln,
  rwm_ln = r_pm_ln,
  rwt_ln = r_pt_ln,
  red_ln = r_ed_ln,
  scg_ln = s_cg_ln,
  swm_ln = s_pm_ln,
  swt_ln = s_pt_ln,
  sed_ln = s_ed_ln,
  lcg_ln = lam_cg_line,
  lwm_ln = lam_pm_line,
  lwt_ln = lam_pt_line,
  led_ln = lam_ed_line,
  pwt_ln = p_wt_ln,
  ped_ln = p_ed_ln
) %>%
  saveRDS(file = "results//figure_marginal_plot_data.rds")

