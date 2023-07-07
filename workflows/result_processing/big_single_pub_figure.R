#Header=========================================================================
# Kenneth Loonam
# July 2023
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
  filter(parameter %in% c("R_cg", "R_dd", "R_wm", "R_wt")) %>%
  mutate(Covariate = case_when(
    parameter == "R_cg" ~ "Cougar Density",
    parameter == "R_dd" ~ "Density Dependence",
    parameter == "R_wm" ~ "PDSI Lag",
    parameter == "R_wt" ~ "PDSI"
  ))

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
    title = "A - Puma density") +
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
    title = "B - Elk density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * 1.1)) +
  xlim(1, 4.8)
