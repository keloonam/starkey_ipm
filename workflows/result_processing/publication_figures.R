#Header=========================================================================
# Kenneth Loonam
# January 2023
# Figures to support IPM publication

#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

# Results to load
load("results//ipm_result_05jan2023_R_pdis.Rdata")
load("data//elk_ipm_data_05jan2023.Rdata")

# Derived result summaries
data <- summary(rslt)
q    <- data$quantiles

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

sc_dat <- q[c(grep("survival_ca", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "Calf Survival") %>%
  filter(year != 1988)

lam_dat <- q[c(grep("lambda", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "lambda") %>%
  filter(year != 1988, year != 2021)

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



r_cov_dat <- cov_dat %>%
  pivot_wider(names_from = covariate) %>%
  mutate(PDSI_lag = lag(PDSI), 'Female Density' = lag(`Female Density`)) %>%
  filter(!is.na(PDSI)) %>%
  full_join(r_dat)
#####
"Start here!!! We are remaking the demographics plot as a multipanel that uses
geom_ribbon to nicely visualize parameters and variance on compact plots.
After that, we are on to the residual plots. Please for the love of god write a 
function for the residuals. No more of this copy pasting bullshit, Kenneth. 
You are better than that."
#Figures========================================================================
dem_fig_text_size <- 5.5
dem_fig_grph_size <- .3

dem_abund <- ggplot(
  data = abundance, 
  aes(x = year, y = mean, color = class, shape = class)) +
  geom_line(size = dem_fig_grph_size) +
  geom_point(size = dem_fig_grph_size * 2) +
  geom_linerange(aes(ymax = uci, ymin = lci), size = dem_fig_grph_size) +
  theme_classic() +
  labs(title = "Elk abundance", y = "N individuals", x = "") +
  ylim(NA,700) +
  scale_color_jco() +
  theme(
    legend.title = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))

dem_calf <- dem_dat %>% filter(class %in% c("Calf Survival", "Recruitment")) %>%
ggplot(data = ., aes(x = year, y = mean, color = class, shape = class)) +
  geom_line(size = dem_fig_grph_size) +
  geom_point(size = dem_fig_grph_size * 2) +
  geom_linerange(aes(ymax = uci, ymin = lci), size = dem_fig_grph_size) +
  theme_classic() +
  labs(title = "Calf demographics", y = "Ratio/Probability", x = "") +
  xlim(1988,NA) +
  scale_color_jco() +
  theme(
    legend.title = element_blank(), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(size = dem_fig_text_size))

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

plot_grid(
  dem_abund, dem_calf, dem_cov, 
  align = "v", 
  ncol = 1,
  label_size = 2)
ggsave(
  "demographics_fig.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 300)

#Posteriors=====================================================================

post_fig_text_size <- 5.5
post_fig_grph_size <- 0.3

ggplot(cov_post, aes(x = value, color = Covariate)) +
  geom_density(size = post_fig_grph_size) +
  geom_vline(xintercept = 0, size = post_fig_grph_size, linetype = "dashed") +
  theme_classic() +
  labs(x = "Value", y = "Density", title = "Recruitment covariate posterior distributions") +
  scale_color_jco() +
  # theme(legend.position = "bottom") +
  # guides(color = guide_legend(nrow = 2, byrow = T)) +
  theme(text = element_text(size = post_fig_text_size))
ggsave("covariate_fig.png", width = 3.5, height = 1.5, units = "in", dpi = 300)

ggplot(cg_recruit, aes(x = Recruitment, color = Cougars)) +
  geom_density(size = post_fig_grph_size) +
  theme_classic() +
  labs(x = "Calves/Cow", y = "Density", title = "Mean recruitment by cougar density") +
  scale_color_jco() +
  # theme(legend.position = "bottom") +
  # guides(color = guide_legend(nrow = 2, byrow = T)) +
  theme(text = element_text(size = post_fig_text_size))
ggsave("real_post_figure.png", width = 3.5, height = 1.5, units = "in", dpi = 300)

#Marginal Plots=================================================================

marg_fig_text_size <- 5.5
marg_fig_grph_size <- 0.05

# Cougars
cougar_marg_plot <- ggplot(data = r_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = r_cov_dat,
    aes(x = `Cougar Index`, y = mean, ymin = lci, ymax = uci),
    size = 0.1
    ) +
  theme_classic() +
  labs(
    x = "Cougar density index (year t)", 
    y = "Calves/Female", 
    title = "A - Cougar density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-1.75, 0.9)

# PDI Lag
pdilag_marg_plot <- ggplot(data = r_pm_line, aes(x = val, y = mci)) +
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
    x = "Palmer drought index (t-1)", 
    y = "Calves/Female", 
    title = "B - PDI prior year") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-1.5, 2)

# PDI Lag
pdi_marg_plot <- ggplot(data = r_pt_line, aes(x = val, y = mci)) +
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
    x = "Palmer drought index (year t)", 
    y = "Calves/Female", 
    title = "C - PDI current year") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-1.5, 2)

# Elk Density Lag
ed_marg_plot <- ggplot(data = r_ed_line, aes(x = val, y = mci)) +
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
    x = "Female elk density (t-1)", 
    y = "Calves/Female", 
    title = "D - Elk density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-2, 2)

plot_grid(
  cougar_marg_plot, pdilag_marg_plot, pdi_marg_plot, ed_marg_plot, 
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//recruitment_marginal_plots_fig.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 300)

