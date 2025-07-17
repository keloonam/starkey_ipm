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
preg_ly <- 0.25

xlim_puma <- c(0, 2.25)
xlim_spei <- c(-2.4, 2.1)
xlim_nelk <- c(0.75, 5.5)

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

