# Kenneth Loonam
# February 2023
# Post hoc analyses of covariate effects on lambda and calf survival

#Environment====================================================================
require(dplyr); require(ggplot2); require(ggsci); require(rjags)
require(cowplot); require(tidyr); require(purrr)
load("results//ipm_result_21apr2023_R&S_pdi.Rdata")
load("data//elk_ipm_data_21apr2023.Rdata")

cg_sd    <- 0.8404219
cg_mn    <- 1.403528
pdi_sd   <- 1.784047
pdi_mn   <- -1.238529
ed_sd    <- 0.9238554
ed_mn    <- 2.993542
#Data Prep======================================================================
lambda_df <- rslt %>%
  map(., as_tibble) %>%
  bind_rows() %>%
  select(contains("lambda")) %>%
  select(-1)
ed_df <- rslt %>%
  map(., as_tibble) %>%
  bind_rows() %>%
  select(contains("N_tot")) %>%
  select(-34)
ed <- ipm_data$elk_density[-34]
cd <- ipm_data$cougar_density[-1]
pdi <- ipm_data$palmer_index[-1]
pdi_lag <- ipm_data$palmer_index[-34]
sca_df <- rslt %>%
  map(., as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_ca"))

lam_quant <- lambda_df %>%
  map(quantile, probs = c(0.025, 0.5, 0.975)) %>%
  bind_rows() %>%
  mutate(
    year = 2:34,
    ed = ed,
    cd = cd,
    pdi = pdi,
    pdi_lag = pdi_lag)

lam_mean <- unlist(map(lambda_df, mean))
ed_mean <- unlist(map(ed_df, mean))

#Analysis=======================================================================

### Lambda ###
ed <- as.numeric(scale(ed_mean))
glm1 <- glm(
  lam_mean ~ cd + pdi + pdi_lag + ed, 
  family = gaussian(link = "log")
)
summary(glm1)
glm1_res <- tibble(
  Covariate = c("Cougar density", "PSDI", "PSDI lag", "Elk density"),
  Estimate = glm1$coefficients[2:5],
  std_err = c(0.025076, 0.023162, 0.024546, 0.026398)
) %>%
  mutate(LCI = Estimate - std_err * 1.96) %>%
  mutate(UCI = Estimate + std_err * 1.96)

#Results prep===================================================================
lam_coef <- summary(glm1)$coefficients %>% as_tibble() %>%
  mutate(
    covariate = c("b0", "cougar", "pdi", "pdi_lag", "elk"),
    mci = Estimate
  ) %>%
  mutate(
    lci = Estimate - `Std. Error` * 1.96,
    uci = Estimate + `Std. Error` * 1.96
  ) %>%
  select(covariate, lci, mci, uci)

x <- seq(-3, 3, length.out = 1000)
lam_cg <- tibble(
  val = x,
  lci = exp(lam_coef$mci[1] + lam_coef$lci[2] * x),
  mci = exp(lam_coef$mci[1] + lam_coef$mci[2] * x),
  uci = exp(lam_coef$mci[1] + lam_coef$uci[2] * x)
)

ilink <- family(glm1)$linkinv

x <- seq(-3, 3, length.out = 1000)

# Cougar CI
n_cg_dat <- tibble(
  cd      = x,
  pdi     = 0,
  pdi_lag = 0,
  ed      = 0
)
lam_cg_line <- as_tibble(predict(glm1, n_cg_dat, se.fit = T)[1:2]) %>%
  mutate(lci = ilink(fit - se.fit * 1.96),
         mci = ilink(fit),
         uci = ilink(fit + se.fit * 1.96),
         val = x * cg_sd + cg_mn) %>%
  select(val, lci, mci, uci)

# PDI lag CI
n_pdilag_dat <- tibble(
  cd      = 0,
  pdi     = 0,
  pdi_lag = x,
  ed      = 0
)
lam_pdilag_line <- as_tibble(predict(glm1, n_pdilag_dat, se.fit = T)[1:2]) %>%
  mutate(lci = ilink(fit - se.fit * 1.96),
         mci = ilink(fit),
         uci = ilink(fit + se.fit * 1.96),
         val = x * pdi_sd + pdi_mn) %>%
  select(val, lci, mci, uci)

# PDI CI
n_pdi_dat <- tibble(
  cd      = 0,
  pdi     = x,
  pdi_lag = 0,
  ed      = 0
)
lam_pdi_line <- as_tibble(predict(glm1, n_pdi_dat, se.fit = T)[1:2]) %>%
  mutate(lci = ilink(fit - se.fit * 1.96),
         mci = ilink(fit),
         uci = ilink(fit + se.fit * 1.96),
         val = x * pdi_sd + pdi_mn) %>%
  select(val, lci, mci, uci)

# ED CI
n_ed_dat <- tibble(
  cd      = 0,
  pdi     = 0,
  pdi_lag = 0,
  ed      = x
)
lam_ed_line <- as_tibble(predict(glm1, n_ed_dat, se.fit = T)[1:2]) %>%
  mutate(lci = ilink(fit - se.fit * 1.96),
         mci = ilink(fit),
         uci = ilink(fit + se.fit * 1.96),
         val = x * ed_sd + ed_mn) %>%
  select(val, lci, mci, uci)

#Figures========================================================================

##### rsf style presentation plot
# Lambda rsf style point range plot
# ggplot(glm1_res, aes(x = Covariate, y = Estimate)) +
#   geom_point() +
#   geom_linerange(ymin = glm1_res$LCI, ymax = glm1_res$UCI) +
#   coord_flip() +
#   theme_classic() +
#   ylim(-.2, .2) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   labs(x = "", title = "Effects on lambda")
# ggsave("lambda_covariate_plot.png", width = 5, height = 3, units = "in", dpi = 300)
#####

# Lambda marginal plots
marg_fig_text_size <- 7
marg_fig_grph_size <- 0.03
# Cougars
cougar_marg_plot <- ggplot(data = lam_cg_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = cd * cg_sd + cg_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = 0.1
  ) +
  theme_classic() +
  labs(
    x = bquote('Puma 100 km'^-2), 
    y = "Lambda", 
    title = "A - Puma density") +
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(0, 2.12)

# PDI lag
pdilag_marg_plot <- ggplot(data = lam_pdilag_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = pdi_lag * pdi_sd + pdi_mn, 
        y = `50%`, 
        ymin = `2.5%`, 
        ymax = `97.5%`),
    size = 0.1
  ) +
  theme_classic() +
  labs(
    x = "PDI (t-1)", 
    y = "Lambda", 
    title = "D - PDI prior year") + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-3.7, 2.3)

# PDI
pdi_marg_plot <- ggplot(data = lam_pdi_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = pdi * pdi_sd + pdi_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = 0.1
  ) +
  theme_classic() +
  labs(
    x = "PDI (t)", 
    y = "Lambda", 
    title = "C - Summer PDI") + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(-4.65, 2.3)

# ED
ed_marg_plot <- ggplot(data = lam_ed_line, aes(x = val, y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = "grey70") +
  geom_line() +
  # geom_point(data = r_cov_dat, aes(x = `Cougar Index`, y = mean)) +
  geom_pointrange(
    data = lam_quant,
    aes(x = ed * ed_sd + ed_mn, y = `50%`, ymin = `2.5%`, ymax = `97.5%`),
    size = 0.1
  ) +
  theme_classic() +
  labs(
    x = bquote('Females km'^-2), 
    y = "Lambda", 
    title = "B - Elk Density") + 
  scale_color_jco() +
  theme(text = element_text(size = marg_fig_text_size)) +
  xlim(1, 4.8)

plot_grid(
  cougar_marg_plot, ed_marg_plot, pdi_marg_plot,pdilag_marg_plot,
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//lambda_marginal_plots_fig.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 300)

#residual plots=================================================================

lam_residuals <- tibble(
  cd = cd,
  ed = ed,
  wt = pdi,
  wm = pdi_lag,
  residual = glm1$residuals
)

resid_fig_text_size <- 4
resid_fig_elem_size <- 0.15

resid_lam_cg <- ggplot(data = lam_residuals, aes(x = cd, y = residual)) +
  geom_point(size = resid_fig_elem_size) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  labs(
    x = "Cougar density index", 
    y = "Lambda residual", 
    title = "A - Cougar index residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_lam_ed <- ggplot(data = lam_residuals, aes(x = ed, y = residual)) +
  geom_point(size = resid_fig_elem_size) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  labs(
    x = "Female abundance", 
    y = "Lambda residual", 
    title = "B - Elk density residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_lam_wt <- ggplot(data = lam_residuals, aes(x = wt, y = residual)) +
  geom_point(size = resid_fig_elem_size) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  labs(
    x = "PDI (year t)", 
    y = "Lambda residual", 
    title = "C - PDSI residuals") +
  theme(text = element_text(size = resid_fig_text_size))

resid_lam_wm <- ggplot(data = lam_residuals, aes(x = wm, y = residual)) +
  geom_point(size = resid_fig_elem_size) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  labs(
    x = "PDSI (year - 1)", 
    y = "Lambda residual", 
    title = "D - PDSI (year - 1) residuals") +
  theme(text = element_text(size = resid_fig_text_size))

plot_grid(
  resid_lam_cg, resid_lam_ed, resid_lam_wt, resid_lam_wm, 
  align = "v", 
  ncol = 2,
  label_size = 2)
ggsave(
  "figures//lambda_residuals_plots_fig.png", 
  width = 4, 
  height = 4, 
  units = "in", 
  dpi = 300)

