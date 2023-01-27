require(tidyverse)
load("results//ipm_result_05jan2023_R_pdis.Rdata")
load("data//elk_ipm_data_05jan2023.Rdata")

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

# x <- 1:nrow(ed_df)
# 
# lm_loop <- function(i){
#   lam <- as.vector(log(as.numeric(lambda_df[i,])))
#   elk <- as.vector(scale(as.numeric(ed_df[i,])))
#   tmp <- lm(lam ~ elk + cd + pdi + pdi_lag)
#   out <- tmp$coefficients
#   return(out)
# }
# 
# cfs <- map(.x = x, .f = lm_loop) %>%
#   bind_rows()

lam_mean <- unlist(map(lambda_df, mean))
ed_mean <- unlist(map(ed_df, mean))
m1 <- lm(lam_mean ~ cd + pdi + pdi_lag + scale(ed_mean))
glm1 <- glm(lam_mean ~ cd + pdi + pdi_lag + scale(ed_mean), family = gaussian(link = "log"))
summary(m1)
summary(glm1)
m1_res <- tibble(
  Covariate = c("Cougar density", "PSDI", "PSDI lag", "Elk density"),
  Estimate = m1$coefficients[2:5],
  std_err = c(0.04346, 0.03974, .041, .04199)
) %>%
  mutate(LCI = Estimate - std_err * 1.96) %>%
  mutate(UCI = Estimate + std_err * 1.96)


ggplot(m1_res, aes(x = Covariate, y = Estimate)) +
  geom_point() +
  geom_linerange(ymin = m1_res$LCI, ymax = m1_res$UCI) +
  coord_flip() +
  theme_classic() + 
  ylim(-.2, .2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", title = "Effects on lambda")
ggsave("lambda_covariate_plot.png", width = 5, height = 3, units = "in", dpi = 300)

sca_mean <- unlist(map(sca_df, mean))
m2 <- lm(sca_mean[-1] ~ cd + pdi + pdi_lag + scale(ed_mean))
glm2 <- glm(sca_mean[-1] ~ cd + pdi + pdi_lag + scale(ed_mean), family = gaussian(link = "logit"))
summary(glm2)
# load("results//ipm_result_11oct2022_R_cgspddinteraction.Rdata")
# mcmcplots::mcmcplot(rslt)
