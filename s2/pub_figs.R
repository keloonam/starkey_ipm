#Header=========================================================================
# Kenneth Loonam
# January 2023
# Figures to support IPM publication

#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags); require(readr)
require(lubridate); require(cowplot); require(tidyr); require(purrr)

# Functions
source("s2//plotting_functions.R")

# Results to load
cv_rs_raw <- readRDS("s2//results//ipmrs_27sep2025_spei12.rds") %>%
  map(as_tibble) %>%
  bind_rows()
nl_rs_raw <- readRDS("s2//results//ipmrs_27sep2025_null.rds") %>%
  map(as_tibble) %>%
  bind_rows()
ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")
preg <- readRDS("results//pregnancy_analysis_13oct2024.rds") %>%
  map(as_tibble) %>%
  bind_rows()
preg_rd <- read_csv("data//pregnancy_rates.csv")

#Data Building==================================================================
# Demographic data
n_ca <- pull_demo(nl_rs_raw, pat = "N_c", yrs = 1988:2023)
n_af <- pull_demo(nl_rs_raw, pat = "N_f", yrs = 1988:2023)
n_am <- pull_demo(nl_rs_raw, pat = "N_m", yrs = 1988:2023)
n_tot <- pull_demo(nl_rs_raw, pat = "N_tot", yrs = 1988:2023)
r_dat <- pull_demo(nl_rs_raw, pat = "R\\[", yrs = 1989:2023)
sf_dat <- pull_demo(nl_rs_raw, pat = "survival_af", yrs = 1989:2023)
sm_dat <- pull_demo(nl_rs_raw, pat = "survival_am", yrs = 1989:2023)
sc_dat <- pull_demo(nl_rs_raw, pat = "survival_ca", yrs = 1989:2023)
lam_dat_n <- build_lfn_tib(nl_rs_raw)
lam_dat_e <- build_el_tib(nl_rs_raw)
p_cov_dat_scaled <- preg %>%
  select(grep("a_preg_lac_prime", names(.))) %>%
  map(quantile, probs = c(.025, .5, .975)) %>%
  bind_rows() %>%
  set_names("lci", "mci", "uci") %>%
  mutate(yr = 1989:2021) %>%
  mutate(spei = ipm_data$spei12[2:34]) %>%
  mutate(nelk = ipm_data$nelk[2:34])

# Covariate data
cov_dat <- tibble(
  covariate = c(
    rep("SPEI", 36), 
    rep("Cougar Index", 36), 
    rep("Female Density", 36)),
  value = c(
    ipm_data$spei12, 
    ipm_data$cdens, 
    ipm_data$nelk),
  year = c(1988:2023, 1988:2023, 1988:2023)
)

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
n_af <- n_af %>%
  mutate(mean_density = mci / starkey_area_km2)
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

net_management <- tibble(
  year = 1988:2023,
  n_mov = apply(ipm_data$n_a_mov + ipm_data$n_c_mov, 2, sum),
  class = c(rep("total", 36))
)

# Beta results
cov_summaries <- cv_rs_raw %>%
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
rb_post <- cv_rs_raw %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_B0", "R_dd", "R_wm", "R_wt")) %>%
  mutate(id = sort(rep(1:(nrow(.)/5), 5))) %>%
  pivot_wider(id_cols = id, names_from = parameter)
sb_post <- cv_rs_raw %>%
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
# Recruitment
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
r_cov_dat_scaled <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(SPEI_lag = lag(SPEI), 'Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(SPEI)) %>%
  mutate(yr = year) %>%
  select(-year) %>%
  full_join(r_dat) %>%
  filter(!is.na(mci))
# Survival
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
s_cov_dat_scaled <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(SPEI_lag = lag(SPEI), 'Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(SPEI)) %>%
  mutate(yr = year) %>%
  select(-year) %>%
  full_join(sc_dat) %>%
  filter(!is.na(mci))
# Lambda
# Puma
cg_lam_dt <- sim_lam(
  fdt = cv_rs_raw, 
  cov = ipm_data$cdens, 
  length_cov = 100,
  r_cov_name = "R_cg",
  s_cov_name = "S_cg")
# SPEI t-1
wm_lam_dt <- sim_lam(
  fdt = cv_rs_raw, 
  cov = ipm_data$spei12, 
  length_cov = 100,
  r_cov_name = "R_wm",
  s_cov_name = "S_wm")
# SPEI
wt_lam_dt <- sim_lam(
  fdt = cv_rs_raw, 
  cov = ipm_data$spei12, 
  length_cov = 50,
  r_cov_name = "R_wt",
  s_cov_name = "S_wt")
# Density dependence
dd_lam_dt <- sim_lam(
  fdt = cv_rs_raw, 
  cov = ipm_data$nelk, 
  length_cov = 50,
  r_cov_name = "R_dd",
  s_cov_name = "S_dd")
# Pregnancy
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

#To real scale
r_cov_dat <- r_cov_dat_scaled %>%
  mutate(
    `Cougar Index` = `Cougar Index` * 
      sd(cougar_dat$cougar_density) + 
      mean(cougar_density)) %>%
  mutate(`SPEI` = `SPEI` * sd(spei_dat$spei) + mean(spei_dat$spei)) %>%
  mutate(`SPEI_lag` = `SPEI_lag` * sd(spei_dat$spei) + mean(spei_dat$spei)) %>%
  mutate(`Female Density` = `Female Density` * 
           sd(n_af$mean_density) + 
           mean(n_af$mean_density)) 
s_cov_dat <- s_cov_dat_scaled %>%
  mutate(
    `Cougar Index` = `Cougar Index` * 
      sd(cougar_dat$cougar_density) + 
      mean(cougar_density)) %>%
  mutate(`SPEI` = `SPEI` * sd(spei_dat$spei) + mean(spei_dat$spei)) %>%
  mutate(`SPEI_lag` = `SPEI_lag` * sd(spei_dat$spei) + mean(spei_dat$spei)) %>%
  mutate(`Female Density` = `Female Density` * 
           sd(n_af$mean_density) + 
           mean(n_af$mean_density)) 
lam_cov_dat <- lam_dat_e %>%
  mutate(cd = r_cov_dat$`Cougar Index`) %>%
  mutate(spei_lag = r_cov_dat$SPEI_lag) %>%
  mutate(spei = r_cov_dat$SPEI) %>%
  mutate(nelk = r_cov_dat$`Female Density`)
p_cov_dat <- p_cov_dat_scaled %>%
  mutate(spei = spei * sd(spei_dat$spei) + mean(spei_dat$spei)) %>%
  mutate(nelk = nelk * 
           sd(n_af$mean_density) + 
           mean(n_af$mean_density)) %>%
  filter(!is.na(mci))
r_cg_ln <- r_cg_line %>%
  mutate(val = val * sd(cougar_dat$cougar_density) + mean(cougar_density))
r_pm_ln <- r_wm_line %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
r_pt_ln <- r_wt_line %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
r_ed_ln <- r_dd_line %>%
  mutate(val = val * sd(n_af$mean_density) + mean(n_af$mean_density))
s_cg_ln <- s_cg_line %>%
  mutate(val = val * sd(cougar_dat$cougar_density) + mean(cougar_density))
s_pm_ln <- s_wm_line %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
s_pt_ln <- s_wt_line %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
s_ed_ln <- s_dd_line %>%
  mutate(val = val * sd(n_af$mean_density) + mean(n_af$mean_density))
lam_cg_line <- cg_lam_dt$lam_df %>%
  mutate(val = val * sd(cougar_dat$cougar_density) + mean(cougar_density))
lam_pm_line <- wm_lam_dt$lam_df %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
lam_pt_line <- wt_lam_dt$lam_df %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
lam_ed_line <- dd_lam_dt$lam_df %>%
  mutate(val = val * sd(n_af$mean_density) + mean(n_af$mean_density))
p_wt_ln <- p_wt_line %>%
  mutate(val = val * sd(spei_dat$spei) + mean(spei_dat$spei))
p_ed_ln <- p_dd_line %>%
  mutate(val = val * sd(n_af$mean_density) + mean(n_af$mean_density))

#Demographic Figures============================================================
dem_fig_text_size <- 11.5
dem_fig_grph_size <- .3
dem_c <- c("#009999", "#95003d", "#0c5f00", "#6e006f")

dem_abund <- n_tot %>% 
  ggplot(
    data = ., 
    aes(x = yr, y = mci)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size, 
    fill = dem_c[1],
    alpha = 0.25) +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = dem_c[1]) +
  theme_classic() +
  labs(title = "a - Elk abundance", y = "N individuals", x = "") +
  ylim(NA,800) +
  scale_color_jco() +
  theme(
    legend.title = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))

dem_recruit <- r_dat %>%
  ggplot(data = ., aes(x = yr, y = mci)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = dem_c[2],
    alpha = 0.25) +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = dem_c[2]) +
  theme_classic() +
  labs(title = "b - Recruitment", y = "Calves / female", x = "") +
  xlim(1988,NA) +
  scale_color_jco() +
  geom_vline(xintercept = 1999, linetype = "dashed", size = dem_fig_grph_size) +
  theme(
    legend.title = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))

dem_ca_surv <- sc_dat %>% 
  ggplot(aes(x = yr, y = mci)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = dem_c[3],
    alpha = 0.25) +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = dem_c[3]) +
  theme_classic() +
  labs(title = "d - Calf survival", y = "Probability", x = "Year") +
  xlim(1988,NA) +
  scale_color_jco() +
  geom_vline(xintercept = 1999, linetype = "dashed", size = dem_fig_grph_size) +
  theme(
    legend.title = element_blank(), 
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))

lambda <- ggplot(data = lam_dat_e, aes(x = yr, y = mci)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = dem_c[4],
    alpha = 0.25) +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = dem_c[4]) +
  theme_classic() +
  labs(title = "c - Population growth", y = "Lambda", x = "Year") +
  xlim(1988,NA) +
  scale_color_jco() +
  geom_hline(yintercept = 1, linetype = "dashed", size = dem_fig_grph_size) +
  theme(
    legend.title = element_blank(), 
    # axis.text.x = element_blank(),
    # axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))


plot_grid(
  dem_abund, dem_recruit, lambda, dem_ca_surv, 
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//demographics_fig.png", 
  width = 6.5, 
  height = 4, 
  units = "in", 
  dpi = 600)

#Covariate Figures==============================================================

cov_fig_text_size <- 11.5
cov_fig_grph_size <- .3

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")

elk_density <- ggplot(
  data = n_af,
  aes(x = yr, y = mean_density)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[2]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "a - Elk density", y = bquote('Elk km'^-2), x = "Year") +
  xlim(1988, NA) +
  scale_color_jco()
puma_density <- ggplot(
  data = cougar_dat,
  aes(x = year, y = cougar_density)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[3]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "b - Puma density", y = bquote('Puma 100 km'^-2), x = "Year") +
  xlim(1988, NA) +
  scale_color_jco()
management_plot <- net_management %>% filter(class == "total") %>%
  ggplot(
    data = .,
    aes(x = year, y = n_mov)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[4]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "c - Management influence", y = "Net elk added", x = "Year") +
  xlim(1988, NA) +
  scale_color_jco()
pdi_plot <- ggplot(
  data = spei_dat,
  aes(x = year, y = spei)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[8]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "d - SPEI (12-month)", y = "SPEI", x = "Year") +
  xlim(1988, NA) +
  scale_color_jco()

plot_grid(
  elk_density, puma_density, management_plot, pdi_plot, 
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//covariates_fig.png", 
  width = 6.5, 
  height = 4, 
  units = "in", 
  dpi = 600)

#Marginal Plot Arguments========================================================

marg_fig_grph_size <- .02
marg_fig_text_size <- 7
tit_mult <- .9
point_mult <- .5
label_color <- "#dddddd"
label_mult <- 0.3

puma_lx <- 0.15
spei_lx <- -2.25
nelk_lx <- 0.85

recr_ly <- 0.105
surv_ly <- 0.045
lamb_ly <- 0.59
preg_ly <- 0.25

xlim_puma <- c(0.1, 2.2)
xlim_spei <- c(-2.3, 2)
xlim_nelk <- c(0.8, 5.2)

ylim_recr <- c(0.07, 0.81)
ylim_surv <- c(0, 1)
ylim_lamb <- c(0.55, 1.4)
ylim_preg <- c(0.15, 1)
mpc <- c("#0000cc", "#cb7500", "#137312", "#5e00a4")
mpa <- 0.5
# blue, orange, green, purple
#Recruitment Marginal Plots=====================================================
# Cougars
r_cougar_marg_plot <- ggplot(data = r_cg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Cougar Index`, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text = element_blank()) +
  xlim(xlim_puma) + ylim(ylim_recr) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.96", 
    x = puma_lx, 
    y = recr_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    alpha = 0.8,
    fill = label_color)
# SPEI t-1
r_pdilag_marg_plot <- ggplot(data = r_pm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = SPEI_lag, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_recr) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.99", 
    x = spei_lx, 
    y = recr_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# SPEI t
r_pdi_marg_plot <- ggplot(data = r_pt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = SPEI, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_recr) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.59", 
    x = spei_lx, 
    y = recr_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# Density dependence
r_ed_marg_plot <- ggplot(data = r_ed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Female Density`, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Calves/Female", 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text.x.bottom = element_blank()) +
  xlim(xlim_nelk) + ylim(ylim_recr) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.08", 
    x = nelk_lx, 
    y = recr_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Survival Marginal Plots========================================================
# Puma
s_cougar_marg_plot <- ggplot(data = s_cg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Cougar Index`, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text = element_blank()) +
  xlim(xlim_puma) + ylim(ylim_surv) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.97", 
    x = puma_lx, 
    y = surv_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# SPEI t-1
s_pdilag_marg_plot <- ggplot(data = s_pm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = SPEI_lag, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_surv) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.91", 
    x = spei_lx, 
    y = surv_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# SPEI
s_pdi_marg_plot <- ggplot(data = s_pt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = SPEI, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_surv) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.39", 
    x = -2.25, 
    y = 0.045, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# Density dependence
s_ed_marg_plot <- ggplot(data = s_ed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Female Density`, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Calf survival", 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text.x.bottom = element_blank()) +
  xlim(xlim_nelk) + ylim(ylim_surv) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.85", 
    x = nelk_lx, 
    y = surv_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Lambda marginal plots==========================================================
# Puma
l_cougar_marg_plot <- ggplot(data = lam_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_cov_dat,
    aes(x = cd, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text = element_blank()) +
  xlim(xlim_puma) + ylim(ylim_lamb) +
  annotate(
    geom = "label", 
    label = "P(d) > 0.99", 
    x = puma_lx, 
    y = lamb_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# SPEI t-1
l_pdilag_marg_plot <- ggplot(data = lam_pm_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  geom_pointrange(
    data = lam_cov_dat,
    aes(x = spei_lag,
        y = mci,
        ymin = lci,
        ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text.x.bottom = element_blank(),
        axis.text.y.left = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_lamb) +
  annotate(
    geom = "label", 
    label = "P(d) > 0.99", 
    x = spei_lx, 
    y = lamb_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# SPEI
l_pdi_marg_plot <- ggplot(data = lam_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_cov_dat,
    aes(x = spei, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_lamb) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.63", 
    x = spei_lx, 
    y = lamb_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# Density Dependence
l_ed_marg_plot <- ggplot(data = lam_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_cov_dat,
    aes(x = nelk, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Population Growth", 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text.x.bottom = element_blank()) +
  xlim(xlim_nelk)+ ylim(ylim_lamb) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.44", 
    x = nelk_lx, 
    y = lamb_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Pregnancy Marginal Plots=======================================================
# SPEI
p_pdi_marg_plot <- ggplot(data = p_wt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[1], alpha = mpa) +
  geom_line(color = mpc[1]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = p_cov_dat,
    aes(x = spei, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote(SPEI[t]), 
    y = NULL, 
    title = NULL) + #bquote("d2 (0.97)")) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_spei) + ylim(ylim_preg) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.98", 
    x = spei_lx, 
    y = preg_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# Density Dependence
p_ed_marg_plot <- ggplot(data = p_ed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[1], alpha = mpa) +
  geom_line(color = mpc[1]) +
  geom_pointrange(
    data = p_cov_dat,
    aes(x = nelk, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = "Pregnancy rate", 
    title = NULL, # bquote("d1 (0.70)")
  ) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_nelk) + ylim(ylim_preg) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.45", 
    x = nelk_lx, 
    y = preg_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)
# Pumas
p_cougar_marg_plot <- ggplot(data = s_cg_ln, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = bquote({"Puma 100 km"^-2}[t]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_puma) + ylim(ylim_preg)
# SPEI t-1
p_pdilag_marg_plot <- ggplot(data = s_cg_line, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = bquote(SPEI[t-1]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult),
        axis.text.y.left = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_preg)

#Plot Grid======================================================================
plot_grid(
  r_ed_marg_plot, r_pdi_marg_plot, r_pdilag_marg_plot, r_cougar_marg_plot,
  s_ed_marg_plot, s_pdi_marg_plot, s_pdilag_marg_plot, s_cougar_marg_plot,
  l_ed_marg_plot, l_pdi_marg_plot, l_pdilag_marg_plot, l_cougar_marg_plot,
  p_ed_marg_plot, p_pdi_marg_plot, p_pdilag_marg_plot, p_cougar_marg_plot,
  align = "v",
  ncol = 4,
  label_size = 2)

ggsave(
  "figures//all_marginal_plots.png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 500)

