#Header=========================================================================
# Kenneth Loonam
# January 2026

#Variables======================================================================
tag <- "full_mean"
folder <- "results//"
run_id <- "fbipm_rslt_06feb2026_"

rslt_file <- paste0(folder, run_id, tag, ".rds")
save_file <- paste0("figures//marginal_plots_wean_", tag, ".png")

# ppptdt %>% pull(lci) %>% min(); ppptdt %>% pull(uci) %>% max()
# snptdt %>% pull(lci) %>% min(); snptdt %>% pull(uci) %>% max()
# scptdt %>% pull(lci) %>% min(); scptdt %>% pull(uci) %>% max()
# sfptdt %>% pull(lci) %>% min(); sfptdt %>% pull(uci) %>% max()

ppcolor <- "#0072B2"
sncolor <- "#D55E00"
sccolor <- "#800080"
sfcolor <- "#009E73"
hccolor <- "#000055"
hfcolor <- "#00ff55"
dpcolor <- "#aa2244"

pplim <- c(0.38, 0.91)
snlim <- c(0.11, 1)
sclim <- c(0.19, 0.96)
sflim <- c(0.20, 1)

cglim <- c(0, 2.15)
ddlim <- c(0.8, 7.1)
wtlim <- c(-2.2, 1.9)

#c(t, r, b, l) - travels clockwise from top
xapm <- c(0,   0.1, 0.2, .1)
yapm <- c(0.2, 0,   -0.2,   0.2)
appm <- c(0.2, 0.1, -.2, -0.1)

cglx <- cglim[1] + 0.2 * (cglim[2] - cglim[1])
ddlx <- ddlim[1] + 0.2 * (ddlim[2] - ddlim[1])
wtlx <- wtlim[1] + 0.2 * (wtlim[2] - wtlim[1])

pply <- pplim[1] + 0.05 * (pplim[2] - pplim[1])
snly <- snlim[1] + 0.05 * (snlim[2] - snlim[1])
scly <- sclim[1] + 0.05 * (sclim[2] - sclim[1])
sfly <- sflim[1] + 0.05 * (sflim[2] - sflim[1])

label_alpha <- 0.5

#cath_cursor
#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

# Results to load
rsraw <- readRDS(rslt_file) %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  mutate(cg = pd_logistic, wt = spei_12m, wm = lag(spei_12m), dd = nelk/77.38) %>%
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
make_label <- function(dtx, prm){
  x <- dtx %>% pull(prm)
  y <- mean(x > 0)
  if(y < 0.5){
    yn <- round(1 - y, 2)
    direction <- "-"
  }else{
    yn <- round(y, 2)
    direction <- "+"
  }
  if(yn == 1){out <- "P > 0.99"}else{out <- paste0("P = ", yn)}
  return(out)
}

# Create holders
emtib <- matrix(data = NA, nrow = length(x), ncol = 4)
emtib[,1] <- x
dimnames(emtib)[[2]] <- c("xv", "lci", "mci", "uci")
ppcgln <- ppddln <- ppwmln <- emtib
sncgln <- snddln <- snwtln <- snwmln <- emtib
sccgln <- scddln <- scwtln <- scwmln <- emtib
sfcgln <- sfddln <- sfwtln <- sfwmln <- emtib
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
ppcgpd <- ppdt %>% make_label("ppcg")
ppddpd <- ppdt %>% make_label("ppdd")
ppwmpd <- ppdt %>% make_label("ppwm")
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
  ppcgln[i,2:4] <- quantile(ycg, c(0.05, 0.5, 0.95))
  ppddln[i,2:4] <- quantile(ydd, c(0.05, 0.5, 0.95))
  ppwmln[i,2:4] <- quantile(ywm, c(0.05, 0.5, 0.95))
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
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
sndt <- rsraw %>%
  select(SNB_cg, SNB_dd, SNB_wt, SNB_wm, SN_B0) 
names(sndt) <- c("sncg", "sndd", "snwt", "snwm", "snb0")
sncgpd <- sndt %>% make_label("sncg")
snddpd <- sndt %>% make_label("sndd")
snwtpd <- sndt %>% make_label("snwt")
snwmpd <- sndt %>% make_label("snwm")
for(i in 1:length(x)){
  ycg <- expit(sndt$snb0 + x[i] * sndt$sncg)
  ydd <- expit(sndt$snb0 + x[i] * sndt$sndd)
  ywt <- expit(sndt$snb0 + x[i] * sndt$snwt)
  ywm <- expit(sndt$snb0 + x[i] * sndt$snwm)
  sncgln[i,2:4] <- quantile(ycg, c(0.05, 0.5, 0.95))
  snddln[i,2:4] <- quantile(ydd, c(0.05, 0.5, 0.95))
  snwtln[i,2:4] <- quantile(ywt, c(0.05, 0.5, 0.95))
  snwmln[i,2:4] <- quantile(ywm, c(0.05, 0.5, 0.95))
}
snlndt <- rbind(sncgln, snddln, snwtln, snwmln) %>%
  as_tibble() %>%
  mutate(covariate = c(
    rep("cg", length(x)), rep("dd", length(x)), 
    rep("wt", length(x)), rep("wm", length(x))
  ))
snptdt <- rsraw %>%
  select(grep("SN\\[", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
scdt <- rsraw %>%
  select(SCB_cg, SCB_dd, SCB_wt, SCB_wm, SC_B0) 
names(scdt) <- c("sccg", "scdd", "scwt", "scwm", "scb0")
sccgpd <- scdt %>% make_label("sccg")
scddpd <- scdt %>% make_label("scdd")
scwtpd <- scdt %>% make_label("scwt")
scwmpd <- scdt %>% make_label("scwm")
for(i in 1:length(x)){
  ycg <- expit(scdt$scb0 + x[i] * scdt$sccg)
  ydd <- expit(scdt$scb0 + x[i] * scdt$scdd)
  ywt <- expit(scdt$scb0 + x[i] * scdt$scwt)
  ywm <- expit(scdt$scb0 + x[i] * scdt$scwm)
  sccgln[i,2:4] <- quantile(ycg, c(0.05, 0.5, 0.95))
  scddln[i,2:4] <- quantile(ydd, c(0.05, 0.5, 0.95))
  scwtln[i,2:4] <- quantile(ywt, c(0.05, 0.5, 0.95))
  scwmln[i,2:4] <- quantile(ywm, c(0.05, 0.5, 0.95))
}
sclndt <- rbind(sccgln, scddln, scwtln, scwmln) %>%
  as_tibble() %>%
  mutate(covariate = c(
    rep("cg", length(x)), rep("dd", length(x)), 
    rep("wt", length(x)), rep("wm", length(x))
  ))
scptdt <- rsraw %>%
  select(grep("SC\\[", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
sfdt <- rsraw %>%
  select(SFB_cg, SFB_dd, SFB_wt, SFB_wm, SF_B0) 
names(sfdt) <- c("sfcg", "sfdd", "sfwt", "sfwm", "sfb0")
sfcgpd <- sfdt %>% make_label("sfcg")
sfddpd <- sfdt %>% make_label("sfdd")
sfwtpd <- sfdt %>% make_label("sfwt")
sfwmpd <- sfdt %>% make_label("sfwm")
for(i in 1:length(x)){
  ycg <- expit(sfdt$sfb0 + x[i] * sfdt$sfcg)
  ydd <- expit(sfdt$sfb0 + x[i] * sfdt$sfdd)
  ywt <- expit(sfdt$sfb0 + x[i] * sfdt$sfwt)
  ywm <- expit(sfdt$sfb0 + x[i] * sfdt$sfwm)
  sfcgln[i,2:4] <- quantile(ycg, c(0.05, 0.5, 0.95))
  sfddln[i,2:4] <- quantile(ydd, c(0.05, 0.5, 0.95))
  sfwtln[i,2:4] <- quantile(ywt, c(0.05, 0.5, 0.95))
  sfwmln[i,2:4] <- quantile(ywm, c(0.05, 0.5, 0.95))
}
sflndt <- rbind(sfcgln, sfddln, sfwtln, sfwmln) %>%
  as_tibble() %>%
  mutate(covariate = c(
    rep("cg", length(x)), rep("dd", length(x)), 
    rep("wt", length(x)), rep("wm", length(x))
  ))
sfptdt <- rsraw %>%
  select(grep("SF\\[", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)


hcptdt <- rsraw %>%
  select(grep("H\\[1,1", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
hfptdt <- rsraw %>%
  select(grep("H\\[2,1", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2023, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
dpptdt <- rsraw %>%
  select(grep("Pf\\[", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  mutate(yr = rep(1989:2022, nrow(rsraw))) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95),
    .by = yr
  ) %>%
  left_join(cvdt)
pp_tim_plt <- ggplot(ppptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ppcolor, alpha = 0.25) +
  geom_line(color = ppcolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "a - Pregnancy")
sn_tim_plt <- ggplot(snptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sncolor, alpha = 0.25) +
  geom_line(color = sncolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "b - Pre-weaning Survival")
sc_tim_plt <- ggplot(scptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sccolor, alpha = 0.25) +
  geom_line(color = sccolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "c - Post-weaning Survival")
sf_tim_plt <- ggplot(sfptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sfcolor, alpha = 0.25) +
  geom_line(color = sfcolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "d - Adult Survival")
hc_tim_plt <- ggplot(hcptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = hccolor, alpha = 0.25) +
  geom_line(color = hccolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "e - Yearling Harvest")
hf_tim_plt <- ggplot(hfptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = hfcolor, alpha = 0.25) +
  geom_line(color = hfcolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "f - Adult Harvest")
dp_tim_plt <- ggplot(dpptdt, aes(x = yr, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = dpcolor, alpha = 0.25) +
  geom_line(color = dpcolor) +
  theme_classic() +
  xlab("Year") +
  ylab("Probability") +
  labs(title = "g - Female Detection Probability")

plot_grid(
  pp_tim_plt, sn_tim_plt, sc_tim_plt, sf_tim_plt, hc_tim_plt, hf_tim_plt, 
  dp_tim_plt,
  ncol = 1
)
ggsave(
  "figures//temporal_plots.tif",
  width = 18, 
  height = 27, 
  scale = 1.2, 
  units = "cm", 
  dpi = 600)



rm(list = ls())
tag <- "full_rcns"
folder <- "results//"
run_id <- "fbipm_rslt_06feb2026_"

rslt_file <- paste0(folder, run_id, tag, ".rds")
save_file <- paste0("figures//marginal_plots_wean_", tag, ".png")

# ppptdt %>% pull(lci) %>% min(); ppptdt %>% pull(uci) %>% max()
# snptdt %>% pull(lci) %>% min(); snptdt %>% pull(uci) %>% max()
# scptdt %>% pull(lci) %>% min(); scptdt %>% pull(uci) %>% max()
# sfptdt %>% pull(lci) %>% min(); sfptdt %>% pull(uci) %>% max()

ppcolor <- "#0072B2"
sncolor <- "#D55E00"
sccolor <- "#800080"
sfcolor <- "#009E73"

pplim <- c(0.38, 0.91)
snlim <- c(0.11, 1)
sclim <- c(0.19, 0.96)
sflim <- c(0.20, 1)

cglim <- c(0, 2.15)
ddlim <- c(0.8, 7.1)
wtlim <- c(-2.2, 1.9)

#c(t, r, b, l) - travels clockwise from top
xapm <- c(0,   0.1, 0.2, .1)
yapm <- c(0.2, 0,   -0.2,   0.2)
appm <- c(0.2, 0.1, -.2, -0.1)

cglx <- cglim[1] + 0.2 * (cglim[2] - cglim[1])
ddlx <- ddlim[1] + 0.2 * (ddlim[2] - ddlim[1])
wtlx <- wtlim[1] + 0.2 * (wtlim[2] - wtlim[1])

pply <- pplim[1] + 0.05 * (pplim[2] - pplim[1])
snly <- snlim[1] + 0.05 * (snlim[2] - snlim[1])
scly <- sclim[1] + 0.05 * (sclim[2] - sclim[1])
sfly <- sflim[1] + 0.05 * (sflim[2] - sflim[1])

label_alpha <- 0.5

#cath_cursor
#Environment====================================================================
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)

# Results to load
rsraw <- readRDS(rslt_file) %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()
cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  filter(yr > 1987)
cvdt$pd_odfw_est[is.na(cvdt$pd_odfw_est)] <- min(cvdt$pd_odfw_est, na.rm = T)
ipm_rslt <- readRDS("results/ipm_rslt_04feb2026_null_fb.rds")

nelk <- ipm_rslt %>% map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>%
  select(grep("Ntot", names(.))) %>% as.matrix() %>% colMeans()
ncow <- ipm_rslt %>% map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>%
  select(grep("N_f\\[", names(.))) %>% as.matrix() %>% colMeans()
cvdt$nelk <- nelk
cgnorm <- readRDS("results//fbipm_rslt_06feb2026_full_norm.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>% 
  select(grep("cdrng", names(.))) %>%
  as.matrix() %>% colMeans()
cvdt$pd_norm <- cgnorm
scl <- function(x){
  out <- (x - mean(x)) / sd(x)
  return(out)
}
nft <- rsraw %>%
  select(grep("N_f\\[", names(.))) %>%
  as.matrix()
nfo <- rsraw %>%
  select(grep("N_fo\\[", names(.))) %>%
  as.matrix()

age <- tmp$annual_age[,3:38] %>% as.matrix()
sex <- tmp$sex %>% pull(sex) %>% (function(x)(x == "F")*1)
kal <- tmp$known_alive[,3:38] %>% as.matrix()
inm <- tmp$herd_assignment[,3:38] %>% as.matrix() %>% (function(x)(x == "main")*1)
tmp <- readRDS("data//capture_handling_data.rds")
old_elk <- (age*sex*kal*inm) %>% 
  as_tibble() %>% 
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(!is.na(val) & val > 0) %>%
  mutate(old = val > 14) %>%
  group_by(yr) %>%
  summarise(
    mean_age = mean(val),
    sd_age = sd(val),
    prop_old = mean(old),
    nobs = n()
  ) %>%
  # print(n = 35) %>%
  ungroup() %>%
  filter(nobs > 50) %>%
  pull(prop_old) 
cvdt$pold <- c(old_elk, NA, NA, NA, NA)
plot_corr <- function(dt, x, y){
  xl <- case_when(
    x == "spei_12m" ~ "SPEI",
    x == "nelk" ~ "Elk Abundance",
    x == "pold" ~ "Proportion Old"
  )
  yl <- case_when(
    y == "pd_reconstruction" ~ "Reconstruction",
    y == "pd_odfw_est" ~ "ODFW Estimate",
    y == "pd_full_mean" ~ "Index Mean",
    y == "pd_logistic" ~ "Logistic Growth",
    y == "spei_12m" ~ "SPEI",
    y == "nelk" ~ "Elk Density",
    y == "pd_norm" ~ "Normal Draw"
  )
  xd <- dt %>% pull(x) %>% scl()
  yd <- dt %>% pull(y) %>% scl()
  r2 <- round(summary(lm(yd ~ xd))$r.squared, 3)
  if(xl == "SPEI"){
    out <- dt %>%
      ggplot(aes(x = get(x), y = get(y))) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() + xlab(xl) + ylab(yl) +
      labs(title = bquote({R^2} == .(r2)))
  }else{
    out <- dt %>%
      ggplot(aes(x = get(x), y = scl(get(y)))) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() + xlab(xl) + ylab(yl) +
      labs(title = bquote({R^2} == .(r2)))
  }
  return(out)
}

spei_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "spei_12m") + xlab(NULL)
nelk_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "nelk") + xlab(NULL)
lgst_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "pd_logistic") + xlab(NULL)
mnes_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "pd_full_mean") + xlab(NULL)
norm_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "pd_norm") + xlab(NULL)
odfw_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "pd_odfw_est") + xlab(NULL)
rcns_cp <- cvdt %>% filter(!is.na(pold)) %>% filter(pold < 0.3) %>% 
  plot_corr(x = "pold", y = "pd_reconstruction")

plot_grid(
  spei_cp, nelk_cp, lgst_cp, mnes_cp, norm_cp, odfw_cp, rcns_cp,
  ncol = 1
)
ggsave("figures//old_elk_corr.png",
       width = 8.5,
       height = 27,
       units = "cm",
       dpi = 600)
