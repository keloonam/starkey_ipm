# Kenneth Loonam
# February 2023
# Post hoc analyses of covariate effects on lambda and calf survival

#Environment====================================================================
require(tidyverse)
load("results//ipm_result_28feb2023_R_pdsi.Rdata")
load("data//elk_ipm_data_05jan2023.Rdata")

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
cd <- ipm_data$cougar_density[-1]
pdi <- ipm_data$palmer_index[-1]
pdi_lag <- ipm_data$palmer_index[-34]
sca_df <- rslt %>%
  map(., as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_ca"))

lam_mean <- unlist(map(lambda_df, mean))
ed_mean <- unlist(map(ed_df, mean))

#Analysis=======================================================================

### Lambda ###
glm1 <- glm(
  lam_mean ~ cd + pdi + pdi_lag + scale(ed_mean), 
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

### Calf Survival ###
sca_mean <- unlist(map(sca_df, mean))
m2 <- lm(sca_mean[-1] ~ cd + pdi + pdi_lag + scale(ed_mean))
glm2 <- glm(sca_mean[-1] ~ cd + pdi + pdi_lag + scale(ed_mean), family = gaussian(link = "logit"))
summary(glm2)

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
#####
"STARTING FROM THIS WEBSITE!!!
  https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/"
#####

#Figures========================================================================

# Lambda rsf style point range plot
ggplot(glm1_res, aes(x = Covariate, y = Estimate)) +
  geom_point() +
  geom_linerange(ymin = glm1_res$LCI, ymax = glm1_res$UCI) +
  coord_flip() +
  theme_classic() + 
  ylim(-.2, .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", title = "Effects on lambda")
ggsave("lambda_covariate_plot.png", width = 5, height = 3, units = "in", dpi = 300)

# Lambda marginal plots
# Cougars
cougar_marg_plot <- ggplot(data = lam_cg, aes(x = val, y = mci)) +
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

#Calf Survival==================================================================


