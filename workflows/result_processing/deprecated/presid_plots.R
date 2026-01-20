# Make the partial residuals / marginal plots figure for the publication
# Updated version as of July 2025
# Kenneth Loonam

# Environment ==================================================================
source("s2//plotting_functions.R")
# Packages
require(dplyr); require(ggplot2); require(ggsci); require(rjags); require(readr)
require(lubridate); require(cowplot); require(tidyr); require(purrr)

# Load the data/results ========================================================
lndt <- readRDS("results//figure_marginal_plot_data.rds")
ptdt <- readRDS("results//figure_partial_residuals.rds")

#Marginal Plot Arguments========================================================

marg_fig_grph_size <- .02
marg_fig_text_size <- 9
tit_mult <- .9
point_mult <- .5
label_color <- "#dddddd"
label_mult <- 0.3

puma_lx <- 0.05
spei_lx <- -2.35
nelk_lx <- 0.8

recr_ly <- 0.105
surv_ly <- 0.045
lamb_ly <- 0.59
preg_ly <- 0.2

xlim_puma <- c(0, 2.25)
xlim_spei <- c(-2.4, 2.1)
xlim_nelk <- c(0.75, 5)

ylim_recr <- c(0.07, 0.81)
ylim_surv <- c(0, 1)
ylim_lamb <- c(0.55, 1.25)
ylim_preg <- c(0.15, 1)
mpc <- c("#0000bb", "#cb7500", "#137312", "#8e00a4")
mpa <- 0.5

pmv <- 0.1

# blue, orange, green, purple
#Recruitment Marginal Plots=====================================================
# Cougars
r_cougar_mplt <- ggplot(data = lndt$rcg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$recruitment,
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
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
    fill = label_color) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI t-1
r_pdilag_mplt <- ggplot(data = lndt$rwm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$recruitment,
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI t
r_pdi_mplt <- ggplot(data = lndt$rwt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$recruitment,
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# Density dependence
r_ed_mplt <- ggplot(data = lndt$red_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[4], alpha = mpa) +
  geom_line(color = mpc[4]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$recruitment,
    aes(x = ed, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))

#Survival Marginal Plots========================================================
# Puma
s_cougar_mplt <- ggplot(data = lndt$scg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$survival,
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI t-1
s_pdilag_mplt <- ggplot(data = lndt$swm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$survival,
    aes(x = wm, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI
s_pdi_mplt <- ggplot(data = lndt$swt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$survival,
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# Density dependence
s_ed_mplt <- ggplot(data = lndt$sed_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[3], alpha = mpa) +
  geom_line(color = mpc[3]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$survival,
    aes(x = ed, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))

#Lambda marginal plots==========================================================
# Puma
l_cougar_mplt <- ggplot(data = lndt$lcg_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$lambda,
    aes(x = cg, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI t-1
l_pdilag_mplt <- ggplot(data = lndt$lwm_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  geom_pointrange(
    data = ptdt$lambda,
    aes(x = wm,
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI
l_pdi_mplt <- ggplot(data = lndt$lwt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$lambda,
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# Density Dependence
l_ed_mplt <- ggplot(data = lndt$led_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[2], alpha = mpa) +
  geom_line(color = mpc[2]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$lambda,
    aes(x = ed, y = mci, ymin = lci, ymax = uci),
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
  xlim(xlim_nelk) + ylim(ylim_lamb) +
  annotate(
    geom = "label", 
    label = "P(d) = 0.44", 
    x = nelk_lx, 
    y = lamb_ly, 
    hjust = 0,
    vjust = 0.5,
    size = marg_fig_text_size * label_mult,
    fill = label_color,
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))

#Pregnancy Marginal Plots=======================================================
# SPEI
p_pdi_mplt <- ggplot(data = lndt$pwt_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[1], alpha = mpa) +
  geom_line(color = mpc[1]) +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ptdt$pregnancy,
    aes(x = wt, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) + #bquote("d2 (0.97)")) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text = element_blank()) +
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# Density Dependence
p_ed_mplt <- ggplot(data = lndt$ped_ln, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = mpc[1], alpha = mpa) +
  geom_line(color = mpc[1]) +
  geom_pointrange(
    data = ptdt$pregnancy,
    aes(x = ed, y = mci, ymin = lci, ymax = uci),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL, # bquote("d1 (0.70)")
  ) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text = element_blank()) +
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
    alpha = 0.8) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# Pumas
p_cougar_mplt <- ggplot(data = lndt$scg_ln, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_puma) + ylim(ylim_preg) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))
# SPEI t-1
p_pdilag_mplt <- ggplot(data = lndt$swm_ln, aes(x = val, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult),
        axis.text = element_blank()) +
  xlim(xlim_spei) + ylim(ylim_preg) +
  theme(plot.margin = unit(c(pmv, pmv, pmv, pmv), "cm"))

#Special axis plots=============================================================
my1 <- .13
my2 <- -0.15
my3 <- .23
my4 <- 0.22

mx1 <- -.1
mx2 <- .1
mx3 <- 0.23
mx4 <- .19

plt_ry <- ggplot(data = lndt$red_ln, aes(x = NULL, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Calves/Female", 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text.x.bottom = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylim(ylim_recr) +
  theme(plot.margin = unit(c(my1, my2, my3, my4), "cm"))
plt_sy <- ggplot(data = lndt$sed_ln, aes(x = NULL, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Calf survival", 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text.x.bottom = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylim(ylim_surv) +
  theme(plot.margin = unit(c(my1, my2, my3, my4), "cm"))
plt_ly <- ggplot(data = lndt$led_ln, aes(x = NULL, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Population growth", 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text.x.bottom = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylim(ylim_lamb) +
  theme(plot.margin = unit(c(my1, my2, my3, my4), "cm"))
plt_py  <- ggplot(data = lndt$ped_ln, aes(x = NULL, y = mci)) +
  theme_classic() +
  labs(
    x = NULL, 
    y = "Pregnancy rate", 
    title = NULL) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_blank(), axis.text.x.bottom = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  ylim(ylim_preg) +
  theme(plot.margin = unit(c(my1, my2, my3, my4), "cm"))

plt_edx <- ggplot(data = lndt$ped_ln, aes(x = val, y = NULL)) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = NULL, 
    title = NULL, # bquote("d1 (0.70)")
  ) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_nelk) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = unit(c(mx1, mx2, mx3, mx4), "cm"))
plt_wtx <- ggplot(data = lndt$pwt_ln, aes(x = val, y = NULL)) +
  theme_classic() +
  labs(
    x = bquote(SPEI[t]), 
    y = NULL, 
    title = NULL) + #bquote("d2 (0.97)")) + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_spei) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = unit(c(mx1, mx2, mx3, mx4), "cm"))
plt_wmx <- ggplot(data = lndt$swm_ln, aes(x = val, y = NULL)) +
  geom_blank() +
  theme_classic() +
  labs(
    x = bquote(SPEI[t-1]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult),
        axis.text.y.left = element_blank()) +
  xlim(xlim_spei) + 
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = unit(c(mx1, mx2, mx3, mx4), "cm"))
plt_cgx <- ggplot(data = lndt$scg_ln, aes(x = val, y = NULL)) +
  theme_classic() +
  labs(
    x = bquote({"Puma 100 km"^-2}[t]), 
    y = NULL, 
    title = NULL) +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size),
        axis.text.y.left = element_blank()) +
  theme(plot.title = element_text(size = marg_fig_text_size * tit_mult)) +
  xlim(xlim_puma) +
  theme(axis.line = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = unit(c(mx1, mx2, mx3, mx4), "cm"))

empty_ <- ggplot() + geom_blank() + theme_classic()
#Plot Grid======================================================================
plot_grid(
  plt_ry, r_ed_mplt, r_pdi_mplt, r_pdilag_mplt, r_cougar_mplt,
  plt_sy, s_ed_mplt, s_pdi_mplt, s_pdilag_mplt, s_cougar_mplt,
  plt_ly, l_ed_mplt, l_pdi_mplt, l_pdilag_mplt, l_cougar_mplt,
  plt_py, p_ed_mplt, p_pdi_mplt, p_pdilag_mplt, p_cougar_mplt,
  empty_, plt_edx,   plt_wtx,    plt_wmx,       plt_cgx,
  rel_widths = c(.25, 1, 1, 1, 1),
  rel_heights = c(1, 1, 1, 1, .25),
  align = "none",
  ncol = 5,
  label_size = 2)

ggsave(
  "figures//all_marginal_plots.png",
  width = 18,
  height = 16,
  units = "cm",
  dpi = 500)

