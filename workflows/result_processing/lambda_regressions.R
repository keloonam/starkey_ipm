require(tidyverse)
load("results//ipm_result_11oct2022_R_null.Rdata")
load("data//elk_ipm_data_26oct2022.Rdata")

lambda_df <- rslt %>%
  map(., as_tibble) %>%
  bind_rows() %>%
  select(contains("lambda"))
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

x <- 1:nrow(ed_df)

lm_loop <- function(i){
  lam <- as.vector(log(as.numeric(lambda_df[i,])))
  elk <- as.vector(scale(as.numeric(ed_df[i,])))
  tmp <- lm(lam ~ elk + cd + pdi + pdi_lag)
  out <- tmp$coefficients
  return(out)
}

cfs <- map(.x = x, .f = lm_loop) %>%
  bind_rows()

lam_mean <- unlist(map(lambda_df, mean))
ed_mean <- unlist(map(ed_df, mean))
m1 <- lm(lam_mean ~ cd + pdi + pdi_lag + scale(ed_mean))
summary(m1)
m1_res <- tibble(
  Covariate = c("Cougar density", "PDI", "PDI lag", "Clk density"),
  Estimate = m1$coefficients[2:5],
  std_err = c(0.0343046, 0.0313621, .0323915, .1670)
) %>%
  mutate(LCI = Estimate - std_err * 1.96) %>%
  mutate(UCI = Estimate + std_err * 1.96)


ggplot(m1_res, aes(x = Covariate, y = Estimate)) +
  geom_linerange(ymin = m1_res$LCI, ymax = m1_res$UCI) +
  coord_flip() +
  theme_classic() + 
  ylim(-.4, .3) +
  geom_hline(yintercept = 0) +
  labs(x = "", title = "Effects on lambda")
ggsave("lambda_covariate_plot.png", width = 5, height = 3, units = "in", dpi = 300)

sca_mean <- unlist(map(sca_df, mean))
m1 <- lm(sca_mean[-1] ~ cd + pdi + pdi_lag + scale(ed_mean))

# load("results//ipm_result_11oct2022_R_cgspddinteraction.Rdata")
# mcmcplots::mcmcplot(rslt)
