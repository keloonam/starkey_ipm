# Environment===================================================================
require(rjags); require(tidyverse)

cntrl <- readRDS("results//submission_2//ipm_comparison_control.rds") %>%
  map(as_tibble) %>%
  bind_rows()
itk <- seq(from = 100, to = nrow(cntrl), length.out = 30000)
cntrl <- cntrl[itk,]
rm(itk)

binom_cts <- readRDS("results//submission_2//binomial_counts.rds") %>%
  map(as_tibble) %>%
  bind_rows()
n_as_cts <- readRDS("results//submission_2//est_n_as_counts.rds") %>%
  map(as_tibble) %>%
  bind_rows()
cts_only <- readRDS("results//submission_2//counts_only.rds") %>%
  map(as_tibble) %>%
  bind_rows()

pull_dem <- function(x, rt_pat, yrs){
  x %>%
    select(grep(rt_pat, names(.))) %>%
    set_names(yrs) %>%
    pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
    group_by(yr) %>%
    summarise(
      lci = quantile(val, 0.025),
      mci = quantile(val, 0.5),
      uci = quantile(val, 0.975),
      var = sd(val)^2
    ) %>%
    ungroup() %>%
    mutate(yr = as.numeric(yr))
}

scl <- function(x){
  ((x - mean(x)) / sd(x)) %>% return()
}

# Demographics==================================================================
cre <- cntrl %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "Control") %>%
  mutate(parameter = "R")
cla <- cntrl %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "Control") %>%
  mutate(parameter = "L")
csc <- cntrl %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "Control") %>%
  mutate(parameter = "SC")
csf <- cntrl %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "Control") %>%
  mutate(parameter = "SF")
cnt <- cntrl %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "Control") %>%
  mutate(parameter = "N")

# Binomial Count Data
bre <- binom_cts %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "Binomial Count") %>%
  mutate(parameter = "R")
bla <- binom_cts %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "Binomial Count") %>%
  mutate(parameter = "L")
bsc <- binom_cts %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "Binomial Count") %>%
  mutate(parameter = "SC")
bsf <- binom_cts %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "Binomial Count") %>%
  mutate(parameter = "SF")
bnt <- binom_cts %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "Binomial Count") %>%
  mutate(parameter = "N")

nacre <- n_as_cts %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "N as Counts") %>%
  mutate(parameter = "R")
nacla <- n_as_cts %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "N as Counts") %>%
  mutate(parameter = "L")
nacsc <- n_as_cts %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "N as Counts") %>%
  mutate(parameter = "SC")
nacsf <- n_as_cts %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "N as Counts") %>%
  mutate(parameter = "SF")
nacnt <- n_as_cts %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "N as Counts") %>%
  mutate(parameter = "N")

core <- cts_only %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "Counts Only") %>%
  mutate(parameter = "R")
cola <- cts_only %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "Counts Only") %>%
  mutate(parameter = "L")
cosc <- cts_only %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "Counts Only") %>%
  mutate(parameter = "SC")
cosf <- cts_only %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "Counts Only") %>%
  mutate(parameter = "SF")
cont <- cts_only %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "Counts Only") %>%
  mutate(parameter = "N")

frs <- bind_rows(cre, csc, csf, cla, cnt,
                 bre, bsc, bsf, bla, bnt,
                 nacre, nacsc, nacsf, nacla, nacnt,
                 core, cosc, cosf, cola, cont)

test_rs <- frs %>%
  filter(model != "Control") 
control_rs <- frs %>%
  filter(model == "Control") %>%
  mutate(c_mci = mci, c_var = var) %>%
  select(yr, parameter, c_mci, c_var)
plot_rs <- test_rs %>%
  select(yr, parameter, lci, mci, uci, var, model) %>%
  full_join(control_rs)

rm <- plot_rs %>%
  filter(model == "Counts Only") %>%
  filter(parameter %in% c("R", "SC", "SF")) %>%
  ggplot(aes(x = c_mci, y = mci, color = parameter, shape = parameter)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab("Recruitment") + labs(title = "Median") +
  theme(legend.title = element_blank(), legend.position = "none")
rs <- plot_rs %>%
  filter(parameter %in% c("R", "SC", "SF")) %>%
  filter(model == "Counts Only") %>%
  ggplot(aes(x = sqrt(c_var), y = sqrt(var), 
             color = parameter, shape = parameter)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) + labs(title = "Standard deviation") +
  theme(legend.title = element_blank(), legend.position = "none")

scm <- plot_rs %>%
  filter(parameter == "SC") %>%
  ggplot(aes(x = c_mci, y = mci, colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab("Calf survival") + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")
scs <- plot_rs %>%
  filter(parameter == "SC") %>%
  ggplot(aes(x = sqrt(c_var), y = sqrt(var), colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

sfm <- plot_rs %>%
  filter(parameter == "SF") %>%
  ggplot(aes(x = c_mci, y = mci, colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab("Female survival") + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")
sfs <- plot_rs %>%
  filter(parameter == "SF") %>%
  ggplot(aes(x = sqrt(c_var), y = sqrt(var), colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

nm <- plot_rs %>%
  filter(parameter == "N") %>%
  ggplot(aes(x = c_mci, y = mci, colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab("Abundance") + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")
ns <- plot_rs %>%
  filter(parameter == "N") %>%
  ggplot(aes(x = sqrt(c_var), y = sqrt(var), colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

lm <- plot_rs %>%
  filter(parameter == "L") %>%
  ggplot(aes(x = c_mci, y = mci, colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab("Lambda") + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")
ls <- plot_rs %>%
  filter(parameter == "L") %>%
  ggplot(aes(x = sqrt(c_var), y = sqrt(var), colour = model, shape = model)) +
  geom_point(position = position_jitter()) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#555555", 
              alpha = 0.5) +
  theme_classic() +
  xlab(NULL) + ylab(NULL) + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "none")

# Betas=========================================================================

pull_betas <- function(x, mdl){
  x %>%
    select(R_CG, R_DD, R_WM, R_WT, SC_CG, SC_DD, SC_WM, SC_WT) %>%
    pivot_longer(cols = 1:ncol(.), names_to = "beta", values_to = "val") %>%
    group_by(beta) %>%
    summarise(
      lci = quantile(val, 0.025),
      mci = quantile(val, 0.5), 
      uci = quantile(val, 0.975),
      var = sd(val)^2) %>%
    ungroup() %>%
    mutate(Model = mdl) %>%
    mutate(Parameter = case_when(
      beta %in% c("R_CG", "R_DD", "R_WM", "R_WT") ~ "Recruitment",
      T ~ "Calf survival"
    )) %>%
    mutate(Covariate = case_when(
      beta %in% c("R_CG", "SC_CG") ~ "Puma density",
      beta %in% c("R_WM", "SC_WM") ~ "SPEI (t-1)",
      beta %in% c("R_WT", "SC_WT") ~ "SPEI (t)",
      T ~ "Elk density"
    )) %>%
    select(Model, Parameter, Covariate, lci, mci, uci, var) %>%
    return()
}
cbd <- cntrl %>% pull_betas(mdl = "Control") %>%
  mutate(c_lci = lci, c_mci = mci, c_uci = uci, c_var = var) %>%
  select(Parameter, Covariate, c_lci, c_mci, c_uci, c_var)
bbd <- binom_cts %>% pull_betas(mdl = "Binomial counts")
nbd <- n_as_cts %>% pull_betas(mdl = "Estimates as counts")
obd <- cts_only %>% pull_betas(mdl = "No estimates")

bd <- bind_rows(bbd, nbd, obd) %>%
  full_join(cbd)

bc_beta_comp <- bd %>%
  filter(Model == "Binomial counts") %>%
  ggplot(aes(x = c_mci, y = mci, color = Parameter)) +
  geom_pointrange(aes(ymax = uci, ymin = lci)) +
  geom_pointrange(aes(xmax = c_uci, xmin = c_lci)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "#444444",
    alpha = 0.5
  ) +
  theme_classic() +
  xlab(NULL) + ylab("Binomial counts") + labs(title = NULL) +
  theme(legend.position = "none")

ne_beta_comp <- bd %>%
  filter(Model == "No estimates") %>%
  ggplot(aes(x = c_mci, y = mci, color = Parameter)) +
  geom_pointrange(aes(ymax = uci, ymin = lci)) +
  geom_pointrange(aes(xmax = c_uci, xmin = c_lci)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "#444444",
    alpha = 0.5
  ) +
  theme_classic() +
  xlab("Control") + ylab("Estimates excluded") + labs(title = NULL) +
  theme(legend.title = element_blank(), legend.position = "bottom")
plot_grid(bc_beta_comp, ne_beta_comp, nrow = 2, ncol = 1, rel_heights = c(2,3))
ggsave("figures//count_type_beta_comp.png", 
       dpi = 600,
       units = "cm",
       height = 12,
       width = 8.5)


trs <- bind_rows(cre, csc, csf, cla, cnt, core, cosc, cosf, cola, cont)

trs %>% filter(parameter == "R") %>%
  ggplot(aes(x = yr, y = mci, color = model)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lci, ymax = uci),
    linetype = 2,
    alpha = 0) +
  xlab("Year") + ylab("Recruitment") +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank())
trs %>% filter(parameter == "SC") %>%
  ggplot(aes(x = yr, y = mci, color = model)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lci, ymax = uci),
    linetype = 2,
    alpha = 0) +
  xlab("Year") + ylab("Calf survival") +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank()) 
trs %>% filter(parameter == "SF") %>%
  ggplot(aes(x = yr, y = mci, color = model)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lci, ymax = uci),
    linetype = 2,
    alpha = 0) +
  xlab("Year") + ylab("Female survival") +
  theme_classic() + 
  theme(legend.position = "none", legend.title = element_blank())
trs %>% filter(parameter == "N") %>%
  ggplot(aes(x = yr, y = mci, color = model)) +
  geom_line() +
  geom_ribbon(
    aes(ymin = lci, ymax = uci),
    linetype = 2,
    alpha = 0) +
  xlab("Year") + ylab("Abundance") +
  theme_classic() + 
  theme(legend.position = "none", legend.title = element_blank())

np <- trs %>%
  mutate(sd = sqrt(var)) %>%
  select(yr, mci, sd, parameter, model) %>%
  pivot_wider(names_from = model, values_from = c(mci, sd)) %>%
  filter(parameter == "N") %>%
  ggplot(aes(x = mci_Control, y = `mci_Counts Only`)) +
  geom_point() +
  theme_classic() +
  xlab("Control") + ylab("Estimates excluded") +
  geom_abline(intercept = 0, slope = 1)
lp <- trs %>%
  mutate(sd = sqrt(var)) %>%
  select(yr, mci, sd, parameter, model) %>%
  pivot_wider(names_from = model, values_from = c(mci, sd)) %>%
  filter(parameter == "L") %>%
  ggplot(aes(x = mci_Control, y = `mci_Counts Only`)) +
  geom_point() +
  theme_classic() +
  xlab("Control") + ylab("Estimates excluded") +
  geom_abline(intercept = 0, slope = 1)



core <- cts_only %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "ee") %>%
  mutate(parameter = "R") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cola <- cts_only %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "ee") %>%
  mutate(parameter = "L") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cosc <- cts_only %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "ee") %>%
  mutate(parameter = "SC") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cosf <- cts_only %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "ee") %>%
  mutate(parameter = "SF") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cont <- cts_only %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "ee") %>%
  mutate(parameter = "N") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)

cre <- cntrl %>%
  pull_dem(rt_pat = "R\\[", yrs = 1989:2023) %>%
  mutate(model = "ctrl") %>%
  mutate(parameter = "R") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cla <- cntrl %>%
  pull_dem(rt_pat = "LAMBDA", yrs = 1989:2023) %>%
  mutate(model = "ctrl") %>%
  mutate(parameter = "L") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
csc <- cntrl %>%
  pull_dem(rt_pat = "SC\\[", yrs = 1989:2023) %>%
  mutate(model = "ctrl") %>%
  mutate(parameter = "SC") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
csf <- cntrl %>%
  pull_dem(rt_pat = "SF\\[", yrs = 1989:2023) %>%
  mutate(model = "ctrl") %>%
  mutate(parameter = "SF") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)
cnt <- cntrl %>%
  pull_dem(rt_pat = "Ntot\\[", yrs = 1988:2023) %>%
  mutate(model = "ctrl") %>%
  mutate(parameter = "N") %>%
  mutate(mci = scl(mci)) %>%
  select(yr, mci, model, parameter)

trs <- bind_rows(cnt, csf, cla, cre, csc, cont, cosf, cola, core, cosc) %>%
  pivot_wider(names_from = model, values_from = mci) 
rp <- trs %>%
  filter(parameter == "R") %>%
  ggplot(aes(x = ctrl, y = ee)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") +
  xlab("Control") + ylab("Estimates Excluded") + labs(title = "Recruitment")
scp <- trs %>%
  filter(parameter == "SC") %>%
  ggplot(aes(x = ctrl, y = ee)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") +
  xlab("Control") + ylab("Estimates Excluded") + labs(title = "Calf Survival")
np <- trs %>%
  filter(parameter == "N") %>%
  ggplot(aes(x = ctrl, y = ee)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") +
  xlab("Control") + ylab("Estimates Excluded") + labs(title = "Abundance")
lp <- trs %>%
  filter(parameter == "L") %>%
  ggplot(aes(x = ctrl, y = ee)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") +
  xlab("Control") + ylab("Estimates Excluded") + labs(title = "Lambda")
sfp <- trs %>%
  filter(parameter == "SF") %>%
  ggplot(aes(x = ctrl, y = ee)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") +
  xlab("Control") + ylab("Estimates Excluded") + labs(title = "Female Survival")

cowplot::plot_grid(rp, scp, sfp, np, lp, nrow = 3, ncol = 2)
ggsave("figures//estimates_excluded_dem_comp.png", dpi = 600, height = 24, width = 18, units = "cm")
