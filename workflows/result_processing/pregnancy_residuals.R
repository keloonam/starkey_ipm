old_lactating_ids <- rd %>%
  filter(age_class == "old") %>%
  filter(lactating == TRUE) %>%
  pull(id)

old_nl_ids <- rd %>%
  filter(age_class == "old") %>%
  filter(lactating == FALSE) %>%
  pull(id)

prime_lactating_ids <- rd %>%
  filter(age_class == "prime") %>%
  filter(lactating == TRUE) %>%
  pull(id)

prime_nl_ids <- rd %>%
  filter(age_class == "prime") %>%
  filter(lactating == FALSE) %>%
  pull(id)

young_ids <- rd %>%
  filter(age_class == "young") %>%
  pull(id)

p_resid <- p_eff_post %>%
  select(grep("rsdl", names(.))) %>%
  pivot_longer(1:ncol(.)) %>%
  group_by(name) %>%
  summarise(
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975),
  ) %>%
  mutate(id = 1:nrow(.)) %>%
  bind_cols(rd) %>%
  mutate(yr = yr - 1988)

plot_data <- tibble(
  yr = 1:length(nimble_constants$pdi),
  pdi = nimble_constants$pdi,
  elk = nimble_constants$elk
) %>% right_join(p_resid) %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(age = age_class) %>%
  mutate(lac = lactating) %>%
  select(id, yr, pdi, elk, lci, mci, uci, age, lac)
  

o_l <- plot_data %>% filter(age == "old") %>% filter(lac == T)
o_n <- plot_data %>% filter(age == "old") %>% filter(lac == F)
p_l <- plot_data %>% filter(age == "prime") %>% filter(lac == T)
p_n <- plot_data %>% filter(age == "prime") %>% filter(lac == F)
y_n <- plot_data %>% filter(age == "young") %>% filter(lac == F)

resid_fig_text_size <- 4
resid_fig_elem_size <- 0.15

o_l_pdi <- ggplot(data = o_l, aes(x = pdi, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "SPEI scaled", 
    y = "Pregnancy residual", 
    title = "Old lactating - SPEI") +
  theme(text = element_text(size = resid_fig_text_size))

o_n_pdi <- ggplot(data = o_n, aes(x = pdi, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "SPEI scaled", 
    y = "Pregnancy residual", 
    title = "Old nonlactating - SPEI") +
  theme(text = element_text(size = resid_fig_text_size))

p_l_pdi <- ggplot(data = p_l, aes(x = pdi, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "SPEI scaled", 
    y = "Pregnancy residual", 
    title = "Prime lactating - SPEI") +
  theme(text = element_text(size = resid_fig_text_size))

p_n_pdi <- ggplot(data = p_n, aes(x = pdi, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "SPEI scaled", 
    y = "Pregnancy residual", 
    title = "Prime nonlactating - SPEI") +
  theme(text = element_text(size = resid_fig_text_size))

y_n_pdi <- ggplot(data = y_n, aes(x = pdi, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "SPEI scaled", 
    y = "Pregnancy residual", 
    title = "Young nonlactating - SPEI") +
  theme(text = element_text(size = resid_fig_text_size))

plot_grid(
  o_l_pdi, o_n_pdi, p_l_pdi, p_n_pdi, y_n_pdi,
  align = "v",
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//preg_pdi_residuals_plots_fig.png",
  width = 4,
  height = 4,
  units = "in",
  dpi = 300)


o_l_elk <- ggplot(data = o_l, aes(x = elk, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Elk density scaled", 
    y = "Pregnancy residual", 
    title = "Old lactating - density") +
  theme(text = element_text(size = resid_fig_text_size))

o_n_elk <- ggplot(data = o_n, aes(x = elk, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Elk density scaled", 
    y = "Pregnancy residual", 
    title = "Old nonlactating - elk") +
  theme(text = element_text(size = resid_fig_text_size))

p_l_elk <- ggplot(data = p_l, aes(x = elk, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Elk density scaled", 
    y = "Pregnancy residual", 
    title = "Prime lactating - elk") +
  theme(text = element_text(size = resid_fig_text_size))

p_n_elk <- ggplot(data = p_n, aes(x = elk, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Elk density scaled", 
    y = "Pregnancy residual", 
    title = "Prime nonlactating - elk") +
  theme(text = element_text(size = resid_fig_text_size))

y_n_elk <- ggplot(data = y_n, aes(x = elk, y = mci)) +
  geom_pointrange(
    aes(ymin = lci, ymax = uci), 
    size = resid_fig_elem_size) +
  geom_hline(yintercept = 0) +
  theme_classic() +
  labs(
    x = "Elk density scaled", 
    y = "Pregnancy residual", 
    title = "Young nonlactating - elk") +
  theme(text = element_text(size = resid_fig_text_size))


plot_grid(
  o_l_elk, o_n_elk, p_l_elk, p_n_elk, y_n_elk,
  align = "v",
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//preg_elk_residuals_plots_fig.png",
  width = 4,
  height = 4,
  units = "in",
  dpi = 300)











############## START HERE ######################################################
# You were in the middle of some quick marginal plots when you realized you have
# to make a bunch of them (high math says 8). So you're going to tackle that 
# tomorrow. It looks like all of the code is in place below, and I just
# confirmed that from here it can plug into the basic pub figure code.
# There are just some formulas to think about, so you have to do that now, sober
# sucker :P



x <- seq(-3, 3, length.out = 1000)


#Prime_dry======================================================================

# PDSI t marginal plot data
pd_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
pd_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
    x[i]*p_eff_post$bpdi_dry_prm #+ 
    # p_eff_post$b_dry_prm_mean
    ) 
  pd_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
pd_pt_line <- as_tibble(pd_pt_line)
names(pd_pt_line) <- c("val", "lci", "mci", "uci")
pd_pt_line$val <- pd_pt_line$val * 1.784047 - 1.238529

# PDSI t marginal plot data
pd_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
pd_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
    x[i]*p_eff_post$bden_dry_prm
    ) 
  pd_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
pd_ed_line <- as_tibble(pd_ed_line)
names(pd_ed_line) <- c("val", "lci", "mci", "uci")
pd_ed_line$val <- pd_ed_line$val * 0.9238554 + 2.993542

#Old_lac========================================================================

# PDSI t marginal plot data
ol_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
ol_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
      x[i]*p_eff_post$bpdi_lac_old + 
      p_eff_post$b_lac_old_mean +
      p_eff_post$b_old_mean
  ) 
  ol_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
ol_pt_line <- as_tibble(ol_pt_line)
names(ol_pt_line) <- c("val", "lci", "mci", "uci")
ol_pt_line$val <- ol_pt_line$val * 1.784047 - 1.238529

# PDSI t marginal plot data
ol_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
ol_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
    x[i]*p_eff_post$bden_lac_old + 
    p_eff_post$b_lac_old_mean +
    p_eff_post$b_old_mean
    ) 
  ol_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
ol_ed_line <- as_tibble(ol_ed_line)
names(ol_ed_line) <- c("val", "lci", "mci", "uci")
ol_ed_line$val <- ol_ed_line$val * 0.9238554 + 2.993542

#Old_dry========================================================================

# PDSI t marginal plot data
od_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
od_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
      x[i]*p_eff_post$bpdi_dry_old + 
      p_eff_post$b_old_mean
  ) 
  od_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
od_pt_line <- as_tibble(od_pt_line)
names(od_pt_line) <- c("val", "lci", "mci", "uci")
od_pt_line$val <- od_pt_line$val * 1.784047 - 1.238529

# PDSI t marginal plot data
od_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
od_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
    x[i]*p_eff_post$bden_dry_old + 
    p_eff_post$b_old_mean) 
  od_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
od_ed_line <- as_tibble(od_ed_line)
names(od_ed_line) <- c("val", "lci", "mci", "uci")
od_ed_line$val <- od_ed_line$val * 0.9238554 + 2.993542

#Young_dry======================================================================

# PDSI t marginal plot data
yd_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
yd_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
      x[i]*p_eff_post$bpdi_dry_yng + 
      p_eff_post$b_yng_mean
  ) 
  yd_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
yd_pt_line <- as_tibble(yd_pt_line)
names(yd_pt_line) <- c("val", "lci", "mci", "uci")
yd_pt_line$val <- yd_pt_line$val * 1.784047 - 1.238529

# PDSI t marginal plot data
yd_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
yd_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(
    p_eff_post$b0_mean + 
      x[i]*p_eff_post$bden_dry_yng + 
      p_eff_post$b_yng_mean) 
  yd_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
yd_ed_line <- as_tibble(yd_ed_line)
names(yd_ed_line) <- c("val", "lci", "mci", "uci")
yd_ed_line$val <- yd_ed_line$val * 0.9238554 + 2.993542

#Point_data_prep================================================================
pdi_mn_now <- r_cov_dat$PDSI %>% mean()
pdi_sd_now <- r_cov_dat$PDSI %>% sd()
elk_mn_now <- n_af %>% pull(mean_density) %>% mean()
elk_sd_now <- n_af %>% pull(mean_density) %>% sd()
point_data_all <- rd %>% 
  mutate(yr = yr - 1988) %>%
  left_join(., cov_dat) %>%
  mutate(pdi = pdi * pdi_sd_now + pdi_mn_now) %>%
  mutate(elk = elk * elk_sd_now + elk_mn_now)

p_all_point_dat <- nimble_results %>%
  map(as_tibble) %>%
  bind_rows() %>%
  map(quantile, probs = c(.025, .5, .975)) %>%
  bind_rows() %>%
  mutate(model_alias = row_names) %>%
  mutate(keep = grepl("p\\[", model_alias)) %>%
  filter(keep == T) %>%
  bind_cols(point_data_all)

pd_pnts <- p_all_point_dat %>%
  filter(age_class == "prime") %>%
  filter(lactating == FALSE)

od_pnts <- p_all_point_dat %>%
  filter(age_class == "old") %>%
  filter(lactating == FALSE)

ol_pnts <- p_all_point_dat %>%
  filter(age_class == "old") %>%
  filter(lactating == TRUE)

yd_pnts <- p_all_point_dat %>%
  filter(age_class == "young") %>%
  filter(lactating == FALSE)

#plots???=======================================================================

marg_fig_grph_size <- .02
marg_fig_text_size <- 7
tit_mult <- .9
point_mult <- .5

pd_ed_marg_plot <- ggplot(data = pd_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = pd_pnts,
    aes(x = elk, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = "Pregnancy rate", 
    title = "Non-lactating, prime-aged", # bquote("d1 (0.70)")
  ) + 
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(0.6, 5.1) + ylim(0, 1) 

od_ed_marg_plot <- ggplot(data = od_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = od_pnts,
    aes(x = elk, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = "Pregnancy rate", 
    title = "Non-lactating, old", # bquote("d1 (0.70)")
  ) + 
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(0.6, 5.1) + ylim(0, 1) 

ol_ed_marg_plot <- ggplot(data = ol_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ol_pnts,
    aes(x = elk, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = "Pregnancy rate", 
    title = "Lactating, old", # bquote("d1 (0.70)")
  ) + 
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(0.6, 5.1) + ylim(0, 1) 

yd_ed_marg_plot <- ggplot(data = yd_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = yd_pnts,
    aes(x = elk, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
  ) +
  theme_classic() +
  labs(
    x = bquote({"Females km"^-2}[t]), 
    y = "Pregnancy rate", 
    title = "Non-lactating, young", # bquote("d1 (0.70)")
  ) + 
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(0.6, 5.1) + ylim(0, 1) 

pd_pdi_marg_plot <- ggplot(data = pd_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = pd_pnts,
    aes(x = pdi, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "SPEI", 
    y = NULL, 
    title = "Non-lactating, prime-aged") +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-2.3, 2) + ylim(0, 1)

od_pdi_marg_plot <- ggplot(data = od_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = od_pnts,
    aes(x = pdi, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "SPEI", 
    y = NULL, 
    title = "Non-lactating, old") +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-2.3, 2) + ylim(0, 1)

ol_pdi_marg_plot <- ggplot(data = ol_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = ol_pnts,
    aes(x = pdi, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "SPEI", 
    y = NULL, 
    title = "Lactating, old") +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-2.3, 2) + ylim(0, 1)

yd_pdi_marg_plot <- ggplot(data = yd_pt_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "#bbbbbb") +
  geom_line(color = "#000000") +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = yd_pnts,
    aes(x = pdi, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = marg_fig_grph_size * point_mult
    # position = position_jitter(width = .2)
  ) +
  theme_classic() +
  labs(
    x = "SPEI", 
    y = NULL, 
    title = "Non-lactating, young") +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-2.3, 2) + ylim(0, 1)

require(cowplot)

plot_grid(
  pd_ed_marg_plot, pd_pdi_marg_plot, od_ed_marg_plot, od_pdi_marg_plot, 
  ol_ed_marg_plot, ol_pdi_marg_plot, yd_ed_marg_plot, yd_pdi_marg_plot,
  align = "v",
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//preg_elk_marginal_plots_fig.png",
  width = 18,
  height = 18,
  units = "cm",
  dpi = 500)
