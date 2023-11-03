#Header=========================================================================
# Kenneth Loonam
# January 2023
# Figures to support IPM publication

#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

# Results to load
load("results//ipm_result_21apr2023_R&S_pdi.Rdata")
load("data//elk_ipm_data_21apr2023.Rdata")
preg <- readRDS("results//pregnancy_analysis_7jul2023.rds")

# Derived result summaries
data <- summary(rslt)
q    <- data$quantiles
which(dimnames(q)[[1]] %in% c("R_mean", "S_C_mean", "S_F_mean", "S_M_mean","LAMBDA_mean"))
q[which(dimnames(q)[[1]] %in% c("R_mean", "S_C_mean", "S_F_mean", "S_M_mean","LAMBDA_mean")),]

#Data Building==================================================================
# Abundance
n_ca <- q[grep("N_c", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Calves")

n_af <- q[grep("N_f", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Adult females")

n_am <- q[grep("N_m", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Adult males")

n_tot <- bind_rows(n_ca, n_af, n_am) %>%
  group_by(year) %>%
  summarise(
    lci = sum(lci),
    mean = sum(mean),
    uci = sum(uci)
  ) %>%
  mutate(class = "Total")

abundance <- bind_rows(n_ca, n_af, n_am, n_tot)

# Recruitment
r_dat <- q[c(grep("R", dimnames(q)[[1]]))[1:33], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1989:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Recruitment")

# Survival
sf_dat <- q[c(grep("survival_af", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Survival - Female")

sm_dat <- q[c(grep("survival_am", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Survival - Male")

sc_dat <- q[c(grep("survival_ca", dimnames(q)[[1]]))[1:33], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1989:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Calf Survival") 

lam_dat <- q[c(grep("lambda", dimnames(q)[[1]]))[1:33], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1989:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "lambda") 

dem_dat <- bind_rows(r_dat, sc_dat, sm_dat, sf_dat)

# Covariates
cov_dat <- tibble(
  covariate = c(
    rep("PDSI", 34), 
    rep("Cougar Index", 34), 
    rep("Female Density", 34)),
  value = c(
    ipm_data$palmer_index, 
    ipm_data$cougar_density, 
    ipm_data$af_density),
  year = c(1988:2021, 1988:2021, 1988:2021)
)

cov_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_dd", "R_wm", "R_wt", "S_cg", "S_dd", "S_wm", "S_wt")) %>%
  mutate(direction = value > 0) %>%
  group_by(parameter) %>%
  summarise(
    lcri = quantile(value, .025),
    medi = quantile(value, 0.50),
    ucri = quantile(value, .975),
    pofd = mean(direction)
  )

# Effects
eff_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_B0", "R_dd", "R_wm", "R_wt")) %>%
  mutate(id = sort(rep(1:(nrow(.)/5), 5))) %>%
  pivot_wider(id_cols = id, names_from = parameter)

expit <- function(x){
  1/(1+exp(-x))
}

#####
high_cg <- expit(eff_post$R_B0 + 0.81*eff_post$R_cg)
low_cg  <- expit(eff_post$R_B0 - 1.6*eff_post$R_cg)
mean_cg <- expit(eff_post$R_B0)
cg_recruit <- tibble(
  Cougars = c(rep("High Density", length(high_cg)), 
              rep("Low Density", length(low_cg)),
              rep("Mean Density", length(mean_cg))
  ),
  Recruitment = c(high_cg, low_cg, mean_cg)
)
#####
x <- seq(-3, 3, length.out = 1000)

# Cougar marginal plot data
r_cg_line <- matrix(data = NA, nrow = length(x), ncol = 4)
r_cg_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(eff_post$R_B0 + x[i]*eff_post$R_cg) 
  r_cg_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
r_cg_line <- as_tibble(r_cg_line)
names(r_cg_line) <- c("val", "lci", "mci", "uci")

# PDSI Lag marginal plot data
r_pm_line <- matrix(data = NA, nrow = length(x), ncol = 4)
r_pm_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(eff_post$R_B0 + x[i]*eff_post$R_wm) 
  r_pm_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
r_pm_line <- as_tibble(r_pm_line)
names(r_pm_line) <- c("val", "lci", "mci", "uci")

# PDSI t marginal plot data
r_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
r_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(eff_post$R_B0 + x[i]*eff_post$R_wt) 
  r_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
r_pt_line <- as_tibble(r_pt_line)
names(r_pt_line) <- c("val", "lci", "mci", "uci")

# PDSI t marginal plot data
r_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
r_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(eff_post$R_B0 + x[i]*eff_post$R_dd) 
  r_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
r_ed_line <- as_tibble(r_ed_line)
names(r_ed_line) <- c("val", "lci", "mci", "uci")



r_cov_dat_scaled <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(PDSI_lag = lag(PDSI), 'Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(PDSI)) %>%
  full_join(r_dat) 

#Demographic Figures============================================================
dem_fig_text_size <- 11.5
dem_fig_grph_size <- .3

dem_abund <- abundance %>% filter(class == "Total") %>% 
  ggplot(
    data = ., 
    aes(x = year, y = mean)) +
    geom_ribbon(
      aes(ymax = uci, ymin = lci), 
      size = dem_fig_grph_size, 
      fill = "#99FFFF") +
    geom_line(size = dem_fig_grph_size, aes(group = 1), colour = "#66CCCC") +
    # geom_point(size = dem_fig_grph_size * 2) +
    theme_classic() +
    labs(title = "a - Elk abundance", y = "N individuals", x = "") +
    ylim(NA,800) +
    scale_color_jco() +
    theme(
      legend.title = element_blank(), 
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      text = element_text(size = dem_fig_text_size))

dem_recruit <- dem_dat %>% filter(class %in% c("Recruitment")) %>%
  ggplot(data = ., aes(x = year, y = mean)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = "#FF9966") +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = "#CC6666") +
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

dem_ca_surv <- dem_dat %>% filter(class %in% c("Calf Survival")) %>%
  ggplot(data = ., aes(x = year, y = mean)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = "#99DDCC") +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = "336600") +
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

lambda <- ggplot(data = lam_dat, aes(x = year, y = mean)) +
  geom_ribbon(
    aes(ymax = uci, ymin = lci), 
    size = dem_fig_grph_size,
    fill = "#FF99CC") +
  geom_line(size = dem_fig_grph_size, aes(group = 1), colour = "#AA99AA") +
  theme_classic() +
  labs(title = "c - Population growth", y = "PP", x = "Year") +
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

dem_cov <- ggplot(
  data = cov_dat, 
  aes(x = year, y = value, color = covariate, shape = covariate)) +
  geom_line(size = dem_fig_grph_size) +
  geom_point(size = dem_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "Recruitment covariates", y = "Centered and Scaled Value", x = "Year") +
  xlim(1988, NA) +
  scale_color_jco() +
  theme(
    legend.title = element_blank(),
    text = element_text(size = dem_fig_text_size))

#Covariate Figures==============================================================

cov_fig_text_size <- 11.5
cov_fig_grph_size <- .3

cougar_density <-  c(0.01488078, 0.02357929, 0.03727105, 0.05868712, 0.09185758,
                     0.14245975, 0.21789379, 0.32656447, 0.47564394, 0.66700050,
                     0.89270325, 1.13381699, 1.36561612, 1.56693693, 1.72693078,
                     1.84534320, 1.92844249, 1.98460714, 2.02160885, 2.04557683,
                     2.06093241, 2.07070110, 2.07688772, 2.08079464, 2.08325745,
                     2.08480819, 2.08578392, 2.08639759, 2.08678343, 2.08702598,
                     2.08717844, 2.08727427, 2.08733449, 2.08737234)
cougar_dat <- tibble(
  year = 1988:2021,
  cougar_density = cougar_density
)
starkey_area_km2 <- 77.382817
n_af <- n_af %>%
  mutate(mean_density = mean / starkey_area_km2)
pdi_real_scale <- c(
  -3.12,  0.04, -1.84, -1.07, -2.44, -0.38, -2.89,  2.19,  1.06,  0.01,  1.88, 
  -2.03, -2.10, -3.49, -3.08, -2.90,  1.12, -1.07, -0.85, -3.56, -0.60, -0.82,
   2.29, -1.16, -1.19,  1.61, -1.26, -2.85, -3.56, -2.00, -2.30, -0.46, -0.76,
  -4.53)
pdi_dat <- tibble(
  year = 1988:2021,
  pdi = pdi_real_scale
)

net_management <- tibble(
  year = rep(1988:2021, 2),
  n_mov = c(apply(ipm_data$n_ad_add - abs(ipm_data$n_ad_rem), 2, sum) +
            apply(ipm_data$n_ca_add - abs(ipm_data$n_ca_rem), 2, sum),
            (ipm_data$n_ad_add - abs(ipm_data$n_ad_rem))[1,] +
            (ipm_data$n_ca_add - abs(ipm_data$n_ca_rem))[1,]),
  class = c(rep("total", 34), rep("females", 34))
)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")

elk_density <- ggplot(
  data = n_af,
  aes(x = year, y = mean_density)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[2]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "a - Elk density", y = bquote('Females km'^-2), x = "Year") +
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
  data = pdi_dat,
  aes(x = year, y = pdi)) +
  geom_line(size = cov_fig_grph_size, color = cbbPalette[8]) +
  geom_point(size = cov_fig_grph_size * 2) +
  theme_classic() +
  labs(title = "d - Palmer index", y = "PDI", x = "Year") +
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

#Marginal Plots=================================================================

marg_fig_text_size <- 7
marg_fig_grph_size <- 0.03

r_cov_dat <- r_cov_dat_scaled %>%
  filter(year != 1988) %>%
  mutate(
    `Cougar Index` = `Cougar Index` * 
                     sd(cougar_dat$cougar_density) + 
                     mean(cougar_density)) %>%
  mutate(`PDSI` = `PDSI` * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)) %>%
  mutate(`PDSI_lag` = `PDSI_lag` * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)) %>%
  mutate(`Female Density` = n_af$mean_density[-34]) 
  
r_cg_ln <- r_cg_line %>%
  mutate(val = val * sd(cougar_dat$cougar_density) + mean(cougar_density))
r_pm_ln <- r_pm_line %>%
  mutate(val = val * sd(pdi_dat$pdi) + mean(pdi_dat$pdi))
r_pt_ln <- r_pt_line %>%
  mutate(val = val * sd(pdi_dat$pdi) + mean(pdi_dat$pdi))
r_ed_ln <- r_ed_line %>%
  mutate(val = val * sd(n_af$mean_density) + mean(n_af$mean_density))

# Cougars
cougar_marg_plot <- ggplot(data = r_cg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Cougar Index`, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    ) +
  theme_classic() +
  labs(
    x = bquote('Puma 100 km'^-2), 
    y = "Calves / female", 
    title = "a - Puma density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(0, 2.12)

# PDI Lag
pdilag_marg_plot <- ggplot(data = r_pm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = PDSI_lag, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "PDI (t-1)", 
    y = NULL, 
    title = bquote(D - PDI[t-1])) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(-3.7, 2.3)

# PDI
pdi_marg_plot <- ggplot(data = r_pt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = PDSI, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "PDI (t)", 
    y = "Calves/Female", 
    title = bquote(C - PDI[t])) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(-4.65, 2.3)

# Elk Density Lag
ed_marg_plot <- ggplot(data = r_ed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Female Density`, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = bquote('Females km'^-2), 
    y = NULL, 
    title = "b - Elk density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(1, 4.8)

# plot_grid(
#   cougar_marg_plot, ed_marg_plot, pdi_marg_plot, pdilag_marg_plot,
#   align = "v", 
#   ncol = 2,
#   label_size = 2)
# ggsave(
#   "figures//recruitment_marginal_plots_fig.png", 
#   width = 4, 
#   height = 4, 
#   units = "in", 
#   dpi = 300)

#Residual plots recruitment=====================================================

r_true_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("R", names(.))) %>%
  select(!grep("R_", names(.))) %>%
  as.matrix()

# eff_post
# r_cov_dat
r_pred_post <- matrix(NA, nrow = nrow(r_true_post), ncol = ncol(r_true_post))
r_residual_post <- matrix(NA, nrow = nrow(r_true_post), ncol = ncol(r_true_post))
for(i in 1:nrow(eff_post)){
  for(j in 2:nrow(r_cov_dat_scaled)){
    r_pred_post[i,j-1] <- expit(
      eff_post$R_B0[i] + 
      eff_post$R_cg[i] *  r_cov_dat_scaled$`Cougar Index`[j] +
      eff_post$R_dd[i] *  r_cov_dat_scaled$`Female Density`[j] +
      eff_post$R_wm[i] *  r_cov_dat_scaled$`PDSI_lag`[j] +
      eff_post$R_wt[i] *  r_cov_dat_scaled$`PDSI`[j])
    r_residual_post[i,j-1] <- r_true_post[i,j-1] - r_pred_post[i,j-1]
  }
}

full_r_dat <- r_residual_post %>%
  as_tibble() %>%
  map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
  bind_rows() %>%
  mutate(
    resid_lci = `2.5%`,
    resid_mci = `50%`,
    resid_uci = `97.5%`
  ) %>%
  select(resid_lci, resid_mci, resid_uci) %>%
  bind_cols(r_cov_dat_scaled[-1,], .)
  
resid_fig_text_size <- 4
resid_fig_elem_size <- 0.15

resid_r_cg <- ggplot(data = full_r_dat, aes(x = `Cougar Index`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Puma density scaled", 
    y = "Recruitment residual", 
    title = "a - Puma density residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_r_ed <- ggplot(data = full_r_dat, aes(x = `Female Density`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Female density scaled", 
    y = "Recruitment residual", 
    title = "b - Density dependence residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_r_wt <- ggplot(data = full_r_dat, aes(x = PDSI, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "PDI scaled (t)", 
    y = "Recruitment residual", 
    title = "c - PDI residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_r_wm <- ggplot(data = full_r_dat, aes(x = `PDSI_lag`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "PDI scaled (t-1)", 
    y = "Recruitment residual", 
    title = "d - PDI (t-1) residuals") +
  theme(text = element_text(size = resid_fig_text_size))

# plot_grid(
#   resid_r_cg, resid_r_ed, resid_r_wt, resid_r_wm, 
#   align = "v", 
#   ncol = 2,
#   label_size = 2)
# ggsave(
#   "figures//recruitment_residuals_plots_fig.png", 
#   width = 4, 
#   height = 4, 
#   units = "in", 
#   dpi = 300)

#Survival_data==================================================================
s_cov_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("S_cg", "S_dd", "S_wm", "S_wt")) %>%
  mutate(Covariate = case_when(
    parameter == "S_cg" ~ "Puma Density",
    parameter == "S_dd" ~ "Density Dependence",
    parameter == "S_wm" ~ "PDI Lag",
    parameter == "S_wt" ~ "PDI"
  ))

logit <- function(x){
  log(x/(1-x))
}
expit <- function(x){
  1/(1+exp(-x))
}

# Effects
s_eff_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("S_cg", "S_C_mean", "S_dd", "S_wm", "S_wt")) %>%
  mutate(id = sort(rep(1:(nrow(.)/5), 5))) %>%
  pivot_wider(id_cols = id, names_from = parameter) %>%
  mutate(S_B0 = logit(S_C_mean))

#####
x <- seq(-3, 3, length.out = 1000)

# Cougar marginal plot data
s_cg_line <- matrix(data = NA, nrow = length(x), ncol = 4)
s_cg_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(s_eff_post$S_B0 + x[i]*s_eff_post$S_cg) 
  s_cg_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
s_cg_line <- as_tibble(s_cg_line)
names(s_cg_line) <- c("val", "lci", "mci", "uci")
s_cg_line$val <- s_cg_line$val * 
                  sd(cougar_dat$cougar_density) + 
                  mean(cougar_dat$cougar_density)

# PDSI Lag marginal plot data
s_pm_line <- matrix(data = NA, nrow = length(x), ncol = 4)
s_pm_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(s_eff_post$S_B0 + x[i]*s_eff_post$S_wm) 
  s_pm_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
s_pm_line <- as_tibble(s_pm_line)
names(s_pm_line) <- c("val", "lci", "mci", "uci")
s_pm_line$val <- s_pm_line$val * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)


# PDSI t marginal plot data
s_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
s_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(s_eff_post$S_B0 + x[i]*s_eff_post$S_wt) 
  s_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
s_pt_line <- as_tibble(s_pt_line)
names(s_pt_line) <- c("val", "lci", "mci", "uci")
s_pt_line$val <- s_pt_line$val * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)

# PDSI t marginal plot data
s_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
s_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(s_eff_post$S_B0 + x[i]*s_eff_post$S_dd) 
  s_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
s_ed_line <- as_tibble(s_ed_line)
names(s_ed_line) <- c("val", "lci", "mci", "uci")
s_ed_line$val <- s_ed_line$val * sd(n_af$mean_density) + mean(n_af$mean_density)

s_cov_dat <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(PDI_lag = lag(PDSI)) %>%
  mutate('Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(PDSI)) %>%
  full_join(sc_dat) %>%
  filter(!is.na(mean)) %>%
  mutate(PDSI = PDSI * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)) %>%
  mutate(PDI_lag = `PDI_lag` * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)) %>%
  mutate(`Female Density` = `Female Density` * 
           sd(n_af$mean_density) + 
           mean(n_af$mean_density)) %>%
  mutate(`Cougar Index` = `Cougar Index` * 
           sd(cougar_dat$cougar_density) + 
           mean(cougar_dat$cougar_density))

s_cov_dat_scaled <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(PDI_lag = lag(PDSI)) %>%
  mutate('Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(PDSI)) %>%
  full_join(sc_dat) %>%
  filter(!is.na(mean)) 

#Survival_marginal_plots========================================================


marg_fig_text_size <- 7
marg_fig_grph_size <- 0.03

# Cougars
s_cougar_marg_plot <- ggplot(data = s_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Cougar Index`, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
  ) +
  theme_classic() +
  labs(
    x = bquote('Puma 100 km'^-2), 
    y = "Calf survival", 
    title = "a - Puma density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(0, 2.12)

# PDI Lag
s_pdilag_marg_plot <- ggplot(data = s_pm_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = PDI_lag, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "PDI (t-1)", 
    y = NULL, 
    title = "d - PDI prior year") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(-3.7, 2.3)

# PDI
s_pdi_marg_plot <- ggplot(data = s_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = PDSI, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "PDI (t)", 
    y = "Calf survival", 
    title = "c - PDI current year") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(-4.65, 2.3)

# Elk Density Lag
s_ed_marg_plot <- ggplot(data = s_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Female Density`, y = mean, ymin = lci, ymax = uci),
    size = marg_fig_grph_size
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = bquote('Females km'^-2), 
    y = NULL, 
    title = "b - Elk density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(1, 4.8)

# plot_grid(
#   s_cougar_marg_plot, s_ed_marg_plot, s_pdi_marg_plot, s_pdilag_marg_plot,
#   align = "v", 
#   ncol = 2,
#   label_size = 2)
# ggsave(
#   "figures//survival_marginal_plots_fig.png", 
#   width = 4, 
#   height = 4, 
#   units = "in", 
#   dpi = 300)

#Residual_plots_survival========================================================
s_true_post <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("survival_ca", names(.))) %>%
  as.matrix()

s_pred_post <- matrix(NA, nrow = nrow(s_true_post), ncol = ncol(s_true_post))
s_residual_post <- matrix(NA, nrow = nrow(s_true_post), ncol = ncol(s_true_post))

for(i in 1:nrow(s_eff_post)){
  for(j in 1:nrow(s_cov_dat_scaled)){
    s_pred_post[i,j] <- expit(
      s_eff_post$S_B0[i] + 
        s_eff_post$S_cg[i] *  s_cov_dat_scaled$`Cougar Index`[j] +
        s_eff_post$S_dd[i] *  s_cov_dat_scaled$`Female Density`[j] +
        s_eff_post$S_wm[i] *  s_cov_dat_scaled$`PDI_lag`[j] +
        s_eff_post$S_wt[i] *  s_cov_dat_scaled$`PDSI`[j])
    s_residual_post[i,j] <- s_true_post[i,j] - s_pred_post[i,j]
  }
}

full_s_dat <- s_residual_post %>%
  as_tibble() %>%
  map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
  bind_rows() %>%
  mutate(
    resid_lci = `2.5%`,
    resid_mci = `50%`,
    resid_uci = `97.5%`
  ) %>%
  select(resid_lci, resid_mci, resid_uci) %>%
  bind_cols(s_cov_dat_scaled, .)

resid_fig_text_size <- 4
resid_fig_elem_size <- 0.15

resid_s_cg <- ggplot(data = full_s_dat, aes(x = `Cougar Index`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Puma density scaled", 
    y = "Recruitment residual", 
    title = "a - Puma density residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_s_ed <- ggplot(data = full_s_dat, aes(x = `Female Density`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Female density scaled", 
    y = "Recruitment residual", 
    title = "b - Density dependence residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_s_wt <- ggplot(data = full_s_dat, aes(x = PDSI, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "PDI scaled (t)", 
    y = "Recruitment residual", 
    title = "c - PDI residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_s_wm <- ggplot(data = full_s_dat, aes(x = `PDI_lag`, y = resid_mci)) +
  geom_pointrange(
    aes(ymin = resid_lci, ymax = resid_uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "PDI scaled (t-1)", 
    y = "Recruitment residual", 
    title = "d - PDI (t-1) residuals") +
  theme(text = element_text(size = resid_fig_text_size))

# plot_grid(
#   resid_s_cg, resid_s_ed, resid_s_wt, resid_s_wm, 
#   align = "v", 
#   ncol = 2,
#   label_size = 2)
# ggsave(
#   "figures//survival_residuals_plots_fig.png", 
#   width = 4, 
#   height = 4, 
#   units = "in", 
#   dpi = 300)

#Plot Assembly==================================================================

marg_fig_grph_size <- .02
marg_fig_text_size <- 7
tit_mult <- .9
point_mult <- .5
label_color <- "#dddddd"
label_mult <- 0.3
##### Recruitment #####
r_color_1 <- "#0033FF"
r_color_2 <- "#00CCFF"

# Cougars
r_cougar_marg_plot <- ggplot(data = r_cg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = r_color_2) +
  geom_line(color = r_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Cougar Index`, y = mean, ymin = lci, ymax = uci),
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
  xlim(0, 2.12) + ylim(0.12, 0.81) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.99", 
    x = .05, 
    y = 0.155, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    alpha = 0.8,
    fill = label_color)

# Recruitment PDI Lag===========================================================
r_pdilag_marg_plot <- ggplot(data = r_pm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = r_color_2) +
  geom_line(color = r_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = PDSI_lag, y = mean, ymin = lci, ymax = uci),
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
  xlim(-3.7, 2.3) + ylim(0.12, 0.81) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.92", 
    x = -3.65, 
    y = 0.155, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

# Recruitment PDI===============================================================
r_pdi_marg_plot <- ggplot(data = r_pt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = r_color_2) +
  geom_line(color = r_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = PDSI, y = mean, ymin = lci, ymax = uci),
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
  xlim(-4.65, 2.3) + ylim(0.12, 0.81) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.58", 
    x = -4.6, 
    y = 0.155, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)


# Recruitment Elk Density Lag===================================================
r_ed_marg_plot <- ggplot(data = r_ed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = r_color_2) +
  geom_line(color = r_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Female Density`, y = mean, ymin = lci, ymax = uci),
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
  xlim(1, 4.8) + ylim(0.12, 0.81) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.44", 
    x = 1.05, 
    y = 0.155, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)


##### Survival #####

s_color_1 <- "#993300"
s_color_2 <- "#FFCC66"
# Survival Cougars==============================================================
s_cougar_marg_plot <- ggplot(data = s_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = s_color_2) +
  geom_line(color = s_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Cougar Index`, y = mean, ymin = lci, ymax = uci),
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
  xlim(0, 2.12) + ylim(0.14, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.99", 
    x = .05, 
    y = 0.185, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Survival PDI Lag===============================================================
s_pdilag_marg_plot <- ggplot(data = s_pm_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = s_color_2) +
  geom_line(color = s_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = PDI_lag, y = mean, ymin = lci, ymax = uci),
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
  xlim(-3.7, 2.3) + ylim(0.14, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.66", 
    x = -3.65, 
    y = 0.185, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Survival PDI===================================================================
s_pdi_marg_plot <- ggplot(data = s_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = s_color_2) +
  geom_line(color = s_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = PDSI, y = mean, ymin = lci, ymax = uci),
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
  xlim(-4.65, 2.3) + ylim(0.14, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.39", 
    x = -4.6, 
    y = 0.185, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

# Survival Elk Density Lag======================================================
s_ed_marg_plot <- ggplot(data = s_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = s_color_2) +
  geom_line(color = s_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = s_cov_dat,
    aes(x = `Female Density`, y = mean, ymin = lci, ymax = uci),
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
  xlim(1, 4.8) + ylim(0.14, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.84", 
    x = 1.05, 
    y = 0.185, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

##### Lambda #####
# You'll need to run analyses from 
# "workflows\results_processing\lambda_regressions_and_plots.R

l_color_1 <- "#006633"
l_color_2 <- "#66CC99"
#Lambda Puma====================================================================
l_cougar_marg_plot <- ggplot(data = lam_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = l_color_2) +
  geom_line(color = l_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = cd * cg_sd + cg_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
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
  xlim(0, 2.12) + ylim(0.5955, 1.4) +
  annotate(
    geom = "label", 
    label = "p = 0.01", 
    x = 0.05, 
    y = 0.635, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Lambda PDI lag=================================================================
l_pdilag_marg_plot <- ggplot(data = lam_pdilag_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = l_color_2) +
  geom_line(color = l_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = pdi_lag * pdi_sd + pdi_mn, 
        y = `50%`, 
        ymin = `2.5%`, 
        ymax = `97.5%`),
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
  xlim(-3.7, 2.3)+ ylim(0.5955, 1.4) +
  annotate(
    geom = "label", 
    label = "p = 0.52", 
    x = -3.65, 
    y = 0.635, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Lambda PDI=====================================================================
l_pdi_marg_plot <- ggplot(data = lam_pdi_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = l_color_2) +
  geom_line(color = l_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = pdi * pdi_sd + pdi_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
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
  xlim(-4.65, 2.3)+ ylim(0.5955, 1.4) +
  annotate(
    geom = "label", 
    label = "p = 0.29", 
    x = -4.6, 
    y = 0.635, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#lambda ED======================================================================
l_ed_marg_plot <- ggplot(data = lam_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = l_color_2) +
  geom_line(color = l_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = ed * ed_sd + ed_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "PP", 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(),
        axis.text.x.bottom = element_blank()) +
  xlim(1, 4.8)+ ylim(0.5955, 1.4) +
  annotate(
    geom = "label", 
    label = "p = 0.18", 
    x = 1.05, 
    y = 0.635, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

##### Pregnancy #####

p_color_1 <- "#660066"
p_color_2 <- "#FF99FF"
#Pregnancy PDI==================================================================
p_pdi_marg_plot <- ggplot(data = p_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = p_color_2) +
  geom_line(color = p_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = p_point_dat,
    aes(x = pdi, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote(PDI[t]), 
    y = NULL, 
    title = NULL) + #bquote("d2 (0.97)")) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(-4.65, 2.3) + ylim(0.25, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.97", 
    x = -4.6, 
    y = 0.30, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)


#Pregnancy Density Dependence===================================================
p_ed_marg_plot <- ggplot(data = p_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = p_color_2) +
  geom_line(color = p_color_1) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = p_point_dat,
    aes(x = elk, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
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
  xlim(1, 4.8) + ylim(0.25, 1) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.70", 
    x = 1.05, 
    y = 0.3, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8)

#Pregnancy Pumas================================================================
p_cougar_marg_plot <- ggplot(data = s_cg_line, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = bquote({"Puma 100 km"^-2}[t]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(0, 2.12) + ylim(0.25, 1)

#Pregnancy PDI Lag==============================================================
p_pdilag_marg_plot <- ggplot(data = s_cg_line, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = bquote(PDI[t-1]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult),
        axis.text.y.left = element_blank()) +
  xlim(-3.7, 2.3) + ylim(0.25, 1)

# save(
#   s_cg_line, 
#   s_ed_line, 
#   s_pt_line, 
#   s_pm_line, 
#   r_cov_dat,
#   p_ed_line,
#   p_point_dat,
#   p_pt_line,
#   lam_ed_line,
#   lam_pdi_line,
#   lam_pdilag_line,
#   lam_cg_line,
#   lam_quant,
#   r_cg_line, 
#   r_ed_line, 
#   r_pt_line, 
#   r_pm_line,
#   s_cov_dat,
#   file = "results//full_plot_data.Rdata"
# )
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

