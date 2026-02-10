#Header=========================================================================
# Kenneth Loonam
# January 2026

#Variables======================================================================
tag <- "full_rcns"
folder <- "results//"
run_id <- "fbipm_rslt_06feb2026_"

rslt_file <- paste0(folder, run_id, tag, ".rds")
save_file <- paste0("figures//marginal_plots_", tag, ".png")

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

#Pregnancy_plots================================================================
ppcgplot <- pplndt %>%
  filter(covariate == "cg") %>%
  mutate(val = unscale(xv, mean(cvdt$cg), sd(cvdt$cg))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ppcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = ppptdt, 
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
    color = ppcolor) +
  theme_classic() +
  xlim(cglim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(pplim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = cglx, y = pply, label = ppcgpd, alpha = label_alpha)

ppddplot <- pplndt %>%
  filter(covariate == "dd") %>%
  mutate(val = unscale(xv, mean(cvdt$dd), sd(cvdt$dd))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ppcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = ppptdt, 
    aes(x = dd, y = mci, ymin = lci, ymax = uci),
    color = ppcolor) +
  theme_classic() +
  xlim(ddlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(pplim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = ddlx, y = pply, label = ppddpd, alpha = label_alpha)

ppwmplot <- pplndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ppcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = ppptdt, 
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
    color = ppcolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(pplim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = pply, label = ppwmpd, alpha = label_alpha)

ppwtplot <- pplndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(pplim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm"))

#Neonate_survival_prep==========================================================
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
  sncgln[i,2:4] <- quantile(ycg, c(0.025, 0.5, 0.975))
  snddln[i,2:4] <- quantile(ydd, c(0.025, 0.5, 0.975))
  snwtln[i,2:4] <- quantile(ywt, c(0.025, 0.5, 0.975))
  snwmln[i,2:4] <- quantile(ywm, c(0.025, 0.5, 0.975))
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
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975),
    .by = yr
  ) %>%
  left_join(cvdt)

#Neonate_survival_plots=========================================================
sncgplot <- snlndt %>%
  filter(covariate == "cg") %>%
  mutate(val = unscale(xv, mean(cvdt$cg), sd(cvdt$cg))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sncolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = snptdt, 
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
    color = sncolor) +
  theme_classic() +
  xlim(cglim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(snlim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = cglx, y = snly, label = sncgpd, alpha = label_alpha)

snddplot <- snlndt %>%
  filter(covariate == "dd") %>%
  mutate(val = unscale(xv, mean(cvdt$dd), sd(cvdt$dd))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sncolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = snptdt, 
    aes(x = dd, y = mci, ymin = lci, ymax = uci),
    color = sncolor) +
  theme_classic() +
  xlim(ddlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(snlim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = ddlx, y = snly, label = snddpd, alpha = label_alpha)

snwtplot <- snlndt %>%
  filter(covariate == "wt") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sncolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = snptdt, 
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
    color = sncolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(snlim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = snly, label = snwtpd, alpha = label_alpha)

snwmplot <- snlndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sncolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = snptdt, 
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
    color = sncolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(snlim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = snly, label = snwmpd, alpha = label_alpha)

#Calf_survival_prep=============================================================
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
  sccgln[i,2:4] <- quantile(ycg, c(0.025, 0.5, 0.975))
  scddln[i,2:4] <- quantile(ydd, c(0.025, 0.5, 0.975))
  scwtln[i,2:4] <- quantile(ywt, c(0.025, 0.5, 0.975))
  scwmln[i,2:4] <- quantile(ywm, c(0.025, 0.5, 0.975))
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
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975),
    .by = yr
  ) %>%
  left_join(cvdt)

#Calf_survival_plots============================================================
sccgplot <- sclndt %>%
  filter(covariate == "cg") %>%
  mutate(val = unscale(xv, mean(cvdt$cg), sd(cvdt$cg))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sccolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = scptdt, 
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
    color = sccolor) +
  theme_classic() +
  xlim(cglim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sclim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = cglx, y = scly, label = sccgpd, alpha = label_alpha)

scddplot <- sclndt %>%
  filter(covariate == "dd") %>%
  mutate(val = unscale(xv, mean(cvdt$dd), sd(cvdt$dd))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sccolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = scptdt, 
    aes(x = dd, y = mci, ymin = lci, ymax = uci),
    color = sccolor) +
  theme_classic() +
  xlim(ddlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sclim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = ddlx, y = scly, label = scddpd, alpha = label_alpha)

scwtplot <- sclndt %>%
  filter(covariate == "wt") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sccolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = scptdt, 
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
    color = sccolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sclim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = scly, label = scwtpd, alpha = label_alpha)

scwmplot <- sclndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sccolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = scptdt, 
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
    color = sccolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sclim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = scly, label = scwmpd, alpha = label_alpha)

#Adult_survival_prep============================================================
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
  sfcgln[i,2:4] <- quantile(ycg, c(0.025, 0.5, 0.975))
  sfddln[i,2:4] <- quantile(ydd, c(0.025, 0.5, 0.975))
  sfwtln[i,2:4] <- quantile(ywt, c(0.025, 0.5, 0.975))
  sfwmln[i,2:4] <- quantile(ywm, c(0.025, 0.5, 0.975))
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
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975),
    .by = yr
  ) %>%
  left_join(cvdt)

#Adult_survival_plots===========================================================
sfcgplot <- sflndt %>%
  filter(covariate == "cg") %>%
  mutate(val = unscale(xv, mean(cvdt$cg), sd(cvdt$cg))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sfcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = sfptdt, 
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
    color = sfcolor) +
  theme_classic() +
  xlim(cglim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sflim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = cglx, y = sfly, label = sfcgpd, alpha = label_alpha)

sfddplot <- sflndt %>%
  filter(covariate == "dd") %>%
  mutate(val = unscale(xv, mean(cvdt$dd), sd(cvdt$dd))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sfcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = sfptdt, 
    aes(x = dd, y = mci, ymin = lci, ymax = uci),
    color = sfcolor) +
  theme_classic() +
  xlim(ddlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sflim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = ddlx, y = sfly, label = sfddpd, alpha = label_alpha)

sfwtplot <- sflndt %>%
  filter(covariate == "wt") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sfcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = sfptdt, 
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
    color = sfcolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sflim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = sfly, label = sfwtpd, alpha = label_alpha)

sfwmplot <- sflndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = sfcolor, alpha = 0.5) +
  geom_line() +
  geom_pointrange(
    data = sfptdt, 
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
    color = sfcolor) +
  theme_classic() +
  xlim(wtlim) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  ylim(sflim) + ylab(element_blank()) + theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(appm, "cm")) +
  annotate("label", x = wtlx, y = sfly, label = sfwmpd, alpha = label_alpha)

#X-axis_plots===================================================================
cgxaplot <- sflndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(cglim) + xlab(bquote("Pumas 100 km"^-2)) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylab(NULL) +
  theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(xapm, "cm"))

ddxaplot <- sflndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(ddlim) + xlab(bquote("Females km"^-2)) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylab(NULL) +
  theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(xapm, "cm"))

wtxaplot <- sflndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(wtlim) + xlab(bquote("SPEI"[t])) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylab(NULL) +
  theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(xapm, "cm"))

wmxaplot <- sflndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(wtlim) + xlab(bquote("SPEI"[t-1])) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylab(NULL) +
  theme(axis.text.y = element_blank()) +
  theme(plot.margin = unit(xapm, "cm"))

#Y-axis_plots===================================================================
ppyaplot <- pplndt %>%
  filter(covariate == "wm") %>%
  mutate(val = unscale(xv, mean(cvdt$wt), sd(cvdt$wt))) %>%
  ggplot(aes(x = val, y = mci)) +
  theme_classic() +
  xlim(c(0,0)) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.ticks = element_blank()) +
  ylim(pplim) + ylab("Pregnancy Probability") +
  theme(axis.line.y = element_blank()) +
  theme(plot.margin = unit(yapm, "cm"))
snyaplot <- snlndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(c(0,0)) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.ticks = element_blank()) +
  ylim(snlim) + ylab("Neonate Survival") +
  theme(axis.line.y = element_blank()) +
  theme(plot.margin = unit(yapm, "cm"))
scyaplot <- sclndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(c(0,0)) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.ticks = element_blank()) +
  ylim(sclim) + ylab("Calf Survival") +
  theme(axis.line.y = element_blank()) +
  theme(plot.margin = unit(yapm, "cm"))
sfyaplot <- sflndt %>%
  ggplot(aes(x = xv, y = mci)) +
  theme_classic() +
  xlim(c(0,0)) + xlab(element_blank()) + theme(axis.text.x = element_blank()) +
  theme(axis.line.x = element_blank(), axis.ticks = element_blank()) +
  ylim(sflim) + ylab("Adult Survival") +
  theme(axis.line.y = element_blank()) +
  theme(plot.margin = unit(yapm, "cm"))
blxyplot <- ggplot() + 
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#Plot_assembly==================================================================

plot_grid(
  ppyaplot, ppddplot, ppwtplot, ppwmplot, ppcgplot,
  snyaplot, snddplot, snwtplot, snwmplot, sncgplot,
  scyaplot, scddplot, scwtplot, scwmplot, sccgplot,
  sfyaplot, sfddplot, sfwtplot, sfwmplot, sfcgplot,
  blxyplot, ddxaplot, wtxaplot, wmxaplot, cgxaplot,
  rel_widths = c(0.25, 1, 1, 1, 1),
  rel_heights = c(1, 1, 1, 1, 0.25),
  align = "none",
  ncol = 5,
  label_size = 2)

ggsave(
  save_file,
  width = 18,
  height = 18,
  scale = 1.2,
  units = "cm",
  dpi = 600)

