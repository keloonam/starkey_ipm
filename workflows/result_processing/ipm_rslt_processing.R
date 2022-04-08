# IPM Results Processing
# Kenneth Loonam
# February 2022

#Environment====================================================================

load("results//ipm_result_25march2022.Rdata")
require(tidyverse); require(rjags)

data <- summary(rslt)

q <- data$quantiles

n_tot <- q[grep("N_tot", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = 1988:2021)
s_af <- q[grep("survival_af", dimnames(q)[[1]]), c(1,3,5)]
s_yf <- q[grep("survival_yf", dimnames(q)[[1]]), c(1,3,5)]
s_am <- q[grep("survival_am", dimnames(q)[[1]]), c(1,3,5)]
s_ym <- q[grep("survival_ym", dimnames(q)[[1]]), c(1,3,5)]
s_ca <- q[grep("survival_ca", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = 1988:2021)
r <- q[grep("R", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = 1988:2021)

ggplot(data = n_tot, aes(x = year, y = `50%`)) +
  geom_line() +
  geom_linerange(aes(ymax = `97.5%`, ymin = `2.5%`)) +
  theme_light() +
  labs(title = "N total", x = "Year", y = "Estimated Elk Population") +
  ylim(0,1000)

ggplot(data = r, aes(x = year, y = `50%`)) +
  geom_line() +
  geom_linerange(aes(ymax = `97.5%`, ymin = `2.5%`)) +
  theme_light() +
  labs(title = "Recruitment", x = "Year", y = "Estimated recruitment rate") +
  ylim(0,1)

ggplot(data = s_ca, aes(x = year, y = `50%`)) +
  geom_line() +
  geom_linerange(aes(ymax = `97.5%`, ymin = `2.5%`)) +
  theme_light() +
  labs(title = "Calf survival", x = "Year", y = "Estimated survival from 0.5 to 1.5 yrs") +
  ylim(0,1)

n_af <- q[grep("N*4,1,", dimnames(q)[[1]]), c(1,3,5)]
ggplot(data = n_af, aes(x = year, y = `50%`)) +
  geom_line() +
  geom_linerange(aes(ymax = `97.5%`, ymin = `2.5%`)) +
  theme_light() +
  labs(title = "Calf survival", x = "Year", y = "Estimated survival from 0.5 to 1.5 yrs") +
  ylim(0,1)
