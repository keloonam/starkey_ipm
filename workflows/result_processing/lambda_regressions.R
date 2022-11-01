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
m1 <- lm(lam_mean ~ cd + pdi + pdi_lag + ed_mean)
summary(m1)

# load("results//ipm_result_11oct2022_R_cgspddinteraction.Rdata")
# mcmcplots::mcmcplot(rslt)
