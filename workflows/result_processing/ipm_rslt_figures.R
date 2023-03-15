# IPM Results Processing
# Kenneth Loonam
# February 2022

#Environment====================================================================

load("results//ipm_result_14mar2023_R_null.Rdata")
require(tidyverse); require(rjags); require(ggsci)

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
  ylim(0,700)
ggsave("elk_abundance_plot.png", width = 5, height = 3, units = "in", dpi = 300)
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

r_dat <- q[c(grep("R", dimnames(q)[[1]]))[1:33], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1989:2021))) %>%
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
ggsave("elk_recruitment_plot.png", width = 5, height = 3, units = "in", dpi = 300)


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
ggsave("calf_survival_plot.png", width = 5, height = 3, units = "in", dpi = 300)

#Posteriors=====================================================================

posteriors <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_dd", "R_wm", "R_wt")) %>%
  mutate(Covariate = case_when(
    parameter == "R_cg" ~ "Cougar Density",
    parameter == "R_dd" ~ "Density Dependence",
    parameter == "R_wm" ~ "Climate Lag",
    parameter == "R_wt" ~ "Climate"
  ))

ggplot(posteriors, aes(x = value, color = Covariate)) +
  geom_density() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  labs(x = "Value", y = "Density", title = "Recruitment Covariate Posteriors - pdsi") +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = T))
ggsave("R_covariate_plot.png", width = 5, height = 3, units = "in", dpi = 300)

posteriors %>%
  filter(Covariate == "Climate Lag") %>%
  summarise(out = mean(value > 0))

#Real_scale_effect==============================================================

posteriors <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(1:ncol(.), names_to = "parameter") %>%
  filter(parameter %in% c("R_cg", "R_B0")) %>%
  mutate(id = sort(rep(1:(nrow(.)/2), 2))) %>%
  pivot_wider(id_cols = id, names_from = parameter)

expit <- function(x){
  1/(1+exp(-x))
}

high_cg <- expit(posteriors$R_B0 + 0.81*posteriors$R_cg)
low_cg <- expit(posteriors$R_B0 - 1.6*posteriors$R_cg)
mean_cg <- expit(posteriors$R_B0)
cg_recruit <- tibble(
  Cougars = c(rep("High Density", length(high_cg)), 
              rep("Low Density", length(low_cg)),
              rep("Mean Density", length(mean_cg))
              ),
  Recruitment = c(high_cg, low_cg, mean_cg)
)
require(ggsci)
ggplot(cg_recruit, aes(x = Recruitment, color = Cougars)) +
  geom_density() +
  theme_classic() +
  labs(x = "Calves/Cow", y = "Density", title = "Recruitment Posteriors by Cougar Density") +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = T))
ggsave("R_cougar_posteriors_plot.png", width = 5, height = 3, units = "in", dpi = 300)

cg_recruit %>%
  group_by(Cougars) %>%
  summarise(r_mean = mean(Recruitment))

#lambda=========================================================================

lam_dat <- q[c(grep("lambda", dimnames(q)[[1]]))[2:34], c(1,3,5)] %>%
  as_tibble() %>%
  mutate(year = sort(rep(1989:2021))) %>%
  group_by(year) %>%
  summarise(
    lci = sum(`2.5%`),
    mean = sum(`50%`),
    uci = sum(`97.5%`)
  ) %>%
  mutate(class = "lambda")

ggplot(data = lam_dat, aes(x = year, y = mean)) +
  geom_line() +
  geom_linerange(aes(ymax = uci, ymin = lci)) +
  theme_light() +
  labs(title = "Calf Survival", y = "Survival Probability", x = "Year")
ggsave("calf_survival_plot.png", width = 5, height = 3, units = "in", dpi = 300)

