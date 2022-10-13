# IPM Results Processing
# Kenneth Loonam
# February 2022

#Environment====================================================================

load("results//ipm_result_11oct2022_R_cgspddinteraction.Rdata")
require(tidyverse); require(rjags)

data <- summary(rslt)

q <- data$quantiles

#Abundance_plots================================================================

n_tot_clean <- q[grep("N_tot", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = 1988:2021)

n_ca <- q[grep("N_c", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "calves")

n_af <- q[grep("N_f", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "adult females")

n_am <- q[grep("N_m", dimnames(q)[[1]]), c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
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

ggplot(data = n_tot, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Total Elk", y = "Abundance", x = "Year") +
  ylim(0,1100)
ggplot(data = n_af, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Adult Female Elk", y = "Abundance", x = "Year") +
  ylim(0,600)
ggplot(data = n_am, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Adult Male Elk", y = "Abundance", x = "Year") +
  ylim(0,500)
ggplot(data = n_ca, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Calves", y = "Abundance", x = "Year") +
  ylim(0,300)

#Recruitment====================================================================

r_dat <- q[c(grep("R", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "recruitment")

ggplot(data = r_dat, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Recruitment", y = "Calves/Cow", x = "Year")

#Survival=======================================================================

sf_dat <- q[c(grep("survival_af", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "female_survival")

ggplot(data = sf_dat, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Female Survival", y = "Survival Probability", x = "Year")

sm_dat <- q[c(grep("survival_am", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "male_survival")

ggplot(data = sm_dat, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Male Survival", y = "Survival Probability", x = "Year")

sc_dat <- q[c(grep("survival_ca", dimnames(q)[[1]]))[1:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1988:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "calf_survival")

ggplot(data = sc_dat, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Calf Survival", y = "Survival Probability", x = "Year")
