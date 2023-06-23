# Pregnancy analysis
# Kenneth Loonam
# May 2023

#Variables======================================================================

niter <- 100000
nburn <- 50000
chain <- 3

#Packages=======================================================================

require(tidyverse); require(nimble)

#Load stuff=====================================================================
rd <- read.csv("data//pregnancy_rates.csv") %>%
  as_tibble()

load("data//elk_ipm_data_21apr2023.Rdata")

source("models//misc//pregnancy.R")

#Prep data======================================================================

cov_dat <- tibble(
  pdi = ipm_data$palmer_index[2:34],
  elk = ipm_data$elk_density[2:34],
  yr  = (1989:2021) - 1988
)

#Prep nimble====================================================================

nimble_data <- list(
  nobs = rd$n_observations,
  preg = round(rd$pregnancy_rate * rd$n_observations)
)  

nimble_constants <- list(
  nr = nrow(rd),
  ny = max(rd$yr) - min(rd$yr) + 1,
  yr = rd$yr - 1988,
  yng  = rd$age_class == "young",
  old  = rd$age_class == "old",
  prm  = rd$age_class == "prime",
  pdi  = cov_dat$pdi,
  elk  = cov_dat$elk,
  lac  = rd$lactating * 1
)

nimble_inits <- function(){
  list(
    b0_mean = rlogis(1, 0, 1),
    b0_sd = runif(1, 0, 10),
    b0 = rlogis(33, 0, 1),
    b_yng_mean = rlogis(1, 0, 1),
    b_old_mean = rlogis(1, 0, 1),
    b_yng_sd = runif(1, 0, 10),
    b_old_sd = runif(1, 0, 10),
    b_yng = rlogis(33, 0, 1),
    b_old = rlogis(33, 0, 1),
    b_den_dry_prm = rlogis(1, 0, 1),
    b_den_dry_yng = rlogis(1, 0, 1),
    b_den_dry_old = rlogis(1, 0, 1),
    b_pdi_dry_prm = rlogis(1, 0, 1),
    b_pdi_dry_yng = rlogis(1, 0, 1),
    b_pdi_dry_old = rlogis(1, 0, 1),
    b_den_lac_prm = rlogis(1, 0, 1),
    b_den_lac_yng = rlogis(1, 0, 1),
    b_den_lac_old = rlogis(1, 0, 1),
    b_pdi_lac_prm = rlogis(1, 0, 1),
    b_pdi_lac_yng = rlogis(1, 0, 1),
    b_pdi_lac_old = rlogis(1, 0, 1)
  )
}

nimble_monitors <- c(
  "b0_mean",
  "b0_sd",
  "b_yng_mean",
  "b_old_mean",
  "b_lac_prm_mean",
  "b_lac_yng_mean",
  "b_lac_old_mean",
  "b_yng_sd",
  "b_old_sd",
  "b_lac_prm_sd",
  "b_lac_yng_sd",
  "b_lac_old_sd",
  "bden_dry_prm",
  "bden_dry_yng",
  "bden_dry_old",
  "bpdi_dry_prm",
  "bpdi_dry_yng",
  "bpdi_dry_old",
  "bden_lac_prm",
  "bden_lac_yng",
  "bden_lac_old",
  "bpdi_lac_prm",
  "bpdi_lac_yng",
  "bpdi_lac_old"
)

#Run model======================================================================

nimble_results <- nimble::nimbleMCMC(
  code = nimble_code,
  constants = nimble_constants,
  data = nimble_data,
  inits = nimble_inits,
  monitors = nimble_monitors,
  niter = niter,
  nburnin = nburn,
  nchains = chain
)
  
mcmcplots::mcmcplot(nimble_results)

#Clean results==================================================================

row_names <- dimnames(nimble_results[[1]])[[2]]

cov_dat <- nimble_results %>%
  map(as_tibble) %>%
  bind_rows() %>%
  map(quantile, probs = c(.025, .5, .975)) %>%
  bind_rows() %>%
  mutate(model_alias = row_names) %>%
  mutate(keep = grepl("bden", model_alias) | grepl("bpdi", model_alias)) %>%
  filter(keep == T) %>%
  mutate(lactating = case_when(
    grepl("dry", model_alias) ~ "Not lactating",
    grepl("lac", model_alias) ~ "Lactating"
  )) %>%
  mutate(age = case_when(
    grepl("old", model_alias) ~ "Old",
    grepl("yng", model_alias) ~ "Young",
    grepl("prm", model_alias) ~ "Prime aged"
  )) %>%
  mutate(covariate = case_when(
    grepl("den", model_alias) ~ "Elk density",
    grepl("pdi", model_alias) ~ "Palmer drought index"
  )) %>%
  mutate(lci = `2.5%`) %>%
  mutate(uci = `97.5%`) %>%
  mutate(med = `50%`) %>%
  select(lactating, age, covariate, lci, med, uci) %>%
  filter(age != "Young" | lactating != "Lactating")

prg_dat <- nimble_results %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(young = b0_mean + b_yng_mean) %>%
  mutate(old_lactating = b0_mean + b_old_mean + b_lac_old_mean) %>%
  mutate(old_open = b0_mean + b_old_mean) %>%
  mutate(prime_lactating = b0_mean + b_lac_prm_mean) %>%
  mutate(prime_open = b0_mean) %>%
  select(prime_open, prime_lactating, old_open, old_lactating, young) %>%
  pivot_longer(cols = everything(), names_to = "Age/lactation class", values_to = "val")  %>%
  mutate(val = exp(val)/(1+exp(val))) %>%
  mutate(`Age/lactation class` = case_when(
    `Age/lactation class` == "old_lactating" ~ "Old - lactating",
    `Age/lactation class` == "old_open" ~ "Old",
    `Age/lactation class` == "prime_lactating" ~ "Prime - lactating",
    `Age/lactation class` == "prime_open" ~ "Prime",
    `Age/lactation class` == "young" ~ "Young",
  ))

#Plot results===================================================================

ggplot(cov_dat, aes(x = covariate, y = med, color = age, shape = lactating)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = .5)) +
  theme(legend.title = element_blank()) +
  labs(x = "", y = "Covariate value", title = "Effects on pregnancy rates") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted")
  
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

ggplot(prg_dat, aes(x = `Age/lactation class`, y = val, color = `Age/lactation class`, fill = `Age/lactation class`)) +
  geom_flat_violin() +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "", y = "Posterior distribution", title = "Mean pregnancy rates") +
  coord_flip()

samp_size <- rd %>%
  mutate(ratio = paste0(
    " ",
    round(pregnancy_rate * n_observations),
    "/",
    n_observations
    )) %>%
  mutate(class = case_when(
    age_class == "old" & lactating == TRUE ~ "Old - lactating",
    age_class == "old" & lactating == FALSE ~ "Old",
    age_class == "prime" & lactating == TRUE ~ "Prime - lactating",
    age_class == "prime" & lactating == FALSE ~ "Prime",
    age_class == "young" & lactating == FALSE ~ "Young",
    age_class == "young" & lactating == TRUE ~ "Young - lactating"
  )) %>%
  select(yr, ratio, class) %>%
  pivot_wider(names_from = class, values_from = ratio) %>%
  write.csv("data//pregnancy_sample_sizes.csv")
