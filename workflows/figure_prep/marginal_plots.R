#Header=========================================================================
# Kenneth Loonam
# January 2026

ppcolor <- "#000000"

#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

# Results to load
rsraw <- readRDS("results//ipm_rslt_05jan2026_all_lgst.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  mutate(cg = pd_logistic, wt = spei_12m, wm = lag(spei_12m), dd = nelk) %>%
  select(yr, dd, wt, wm, cg) %>%
  filter(yr > 1987)

# Vector to use for marginal plot calculations
x <- seq(-2.1, 2.25, length.out = 500)


# Functions
expit <- function(x){
  1/(1+exp(-x))
}
unscale <- function(x, org_mn = NA, org_sd = NA){
  return(x * org_sd + org_mn)
}

# Create holders
emtib <- matrix(data = NA, nrow = length(x), ncol = 4)
emtib[,1] <- x
dimnames(emtib)[[2]] <- c("xv", "lci", "mci", "uci")
ppcgln <- ppddln <- ppwmln <- emtib
#Pregnancy_prep=================================================================
ppdt <- rsraw %>%
  select(PB_cg, PB_dd, PB_wm, `P_B0[1]`, `P_B0[2]`, `P_B0[3]`) 
ppyrdt <- rsraw %>%
  select(grep("P_Byr\\[", names(.)), grep("N_f", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(
    c(2:36, 1:36, 1:36, 1:36, 1:36), nrow(rsraw)
  )) %>%
  mutate(param = rep(
    c(rep("ppyr", 35), rep("nft", 36), rep("nfo", 36), rep("nfp", 36), rep("nfy", 36)), 
    nrow(rsraw)
  )) %>%
  mutate(stpn = sort(rep(1:nrow(rsraw), 179))) %>%
  select(stpn, yr, param, value) %>%
  pivot_wider(names_from = param, values_from = value, id_cols = c(stpn, yr)) %>%
  mutate(
    pryf = nfy / nft,
    prpf = nfp / nft,
    prof = nfo / nft
  ) 
mn_age_dist <- ppyrdt %>%
  summarise(pryf = mean(pryf), prpf = mean(prpf), prof = mean(prof))
names(ppdt) <- c("ppcg", "ppdd", "ppwm", "ppb0y", "ppb0p", "ppb0o")
for(i in 1:length(x)){
  ycg <- expit(ppdt$ppb0y + x[i] * ppdt$ppcg) * mn_age_dist$pryf +
         expit(ppdt$ppb0p + x[i] * ppdt$ppcg) * mn_age_dist$prpf +
         expit(ppdt$ppb0o + x[i] * ppdt$ppcg) * mn_age_dist$prof 
  ydd <- expit(ppdt$ppb0y + x[i] * ppdt$ppdd) * mn_age_dist$pryf +
         expit(ppdt$ppb0p + x[i] * ppdt$ppdd) * mn_age_dist$prpf +
         expit(ppdt$ppb0o + x[i] * ppdt$ppdd) * mn_age_dist$prof 
  ywm <- expit(ppdt$ppb0y + x[i] * ppdt$ppwm) * mn_age_dist$pryf +
         expit(ppdt$ppb0p + x[i] * ppdt$ppwm) * mn_age_dist$prpf +
         expit(ppdt$ppb0o + x[i] * ppdt$ppwm) * mn_age_dist$prof 
  ppcgln[i,2:4] <- quantile(ycg, c(0.025, 0.5, 0.975))
  ppddln[i,2:4] <- quantile(ydd, c(0.025, 0.5, 0.975))
  ppwmln[i,2:4] <- quantile(ywm, c(0.025, 0.5, 0.975))
}
pplndt <- rbind(ppcgln, ppddln, ppwmln) %>%
  as_tibble() %>%
  mutate(covariate = c(
    rep("cg", length(x)), rep("dd", length(x)), rep("wm", length(x))
  ))
ppptdt <- rsraw %>%
  select(grep("Pbar", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975),
    .by = yr
  ) %>%
  left_join(cvdt)


pplndt %>%
  filter(covariate == "cg") %>%
  mutate(val = unscale(xv, mean(cvdt$cg), sd(cvdt$cg))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ppcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = ppptdt, 
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
    color = ppcolor)


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
