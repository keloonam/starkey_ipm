# IPM Results Processing
# Kenneth Loonam
# February 2022

#Environment====================================================================

load("results//ipm_result_12apr2022.Rdata")
require(tidyverse); require(rjags)

data <- summary(rslt)

q <- data$quantiles

n_tot_clean <- q[grep("N_tot", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = 1988:2021)
# s_af <- q[grep("survival_af", dimnames(q)[[1]]), c(1,3,5)]
# s_yf <- q[grep("survival_yf", dimnames(q)[[1]]), c(1,3,5)]
# s_am <- q[grep("survival_am", dimnames(q)[[1]]), c(1,3,5)]
# s_ym <- q[grep("survival_ym", dimnames(q)[[1]]), c(1,3,5)]
# s_ca <- q[grep("survival_ca", dimnames(q)[[1]]), c(1,3,5)] %>%
#   as_tibble() %>%
#   mutate(year = 1988:2021)
# r <- q[grep("R", dimnames(q)[[1]]), c(1,3,5)] 
  
n_ca <- q[grep("N.1,", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021,2))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "calves")

n_af <- q[grep("N.[2,3,4],1,", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021,3))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "adult females")

n_am <- q[grep("N.[2,3,4],2,", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021,3))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "adult males")

n_tot <- bind_rows(n_ca, n_af, n_am) %>%
  group_by(year) %>%
  summarise(
    lci = sum(lci),
    mean = sum(mean),
    uci = sum(uci)
  ) %>%
  mutate(class = "total")

abundance <- bind_rows(n_ca, n_af, n_am, n_tot)

ggplot(data = abundance, aes(x = year, y = mean, color = class)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light()

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
