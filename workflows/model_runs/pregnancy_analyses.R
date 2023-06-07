# Pregnancy analysis
# Kenneth Loonam
# May 2023

#Variables======================================================================

niter <- 50000
nburn <- 10000
chain <- 3

#Packages=======================================================================

require(tidyverse); require(nimble)

#Load stuff=====================================================================
rd <- read.csv("data//pregnancy_lactation_2002_present.csv") %>%
  as_tibble()

load("data//elk_ipm_data_21apr2023.Rdata")

source("models//misc//pregnancy.R")

#Prep data======================================================================

preg_rates <- rd %>% group_by(HandlingYr, Lactating) %>%
  summarise(
    pregnancy_rate = mean(Pregnant, na.rm = T),
    pregnancy_sd = sd(Pregnant, na.rm = T),
    n_observations = n()
    ) %>%
  filter(Lactating == "F" | Lactating == "T") %>%
  mutate(lactating = Lactating == "T") %>%
  ungroup() %>%
  separate(HandlingYr, c("year", "yr2")) %>%
  select(year, lactating, pregnancy_rate, pregnancy_sd, n_observations) %>%
  mutate(id = 1:nrow(.))

full_dat <- tibble(
  pdi = rep(ipm_data$palmer_index[15:34], 2),
  elk = rep(ipm_data$elk_density[15:34], 2),
  yr = rep(2002:2021, 2)
) %>%
  arrange(yr) %>%
  mutate(id = 1:nrow(.)) %>%
  full_join(preg_rates)

#Prep nimble====================================================================

nimble_data <- list(
  lac  = full_dat$lactating * 1,
  nobs = full_dat$n_observations,
  preg = round(full_dat$pregnancy_rate * full_dat$n_observations),
  pdi  = full_dat$pdi,
  elk  = full_dat$elk
)  

nimble_constants <- list(
  nr = nrow(full_dat),
  ny = length(unique(full_dat$year)),
  yr   = sort(rep(1:length(unique(full_dat$year)), 2))
)

nimble_inits <- function(){
  list(
    b0_mean = rlogis(1, 0, 1),
    b0_sd = runif(1, 0, 10),
    b0 = rlogis(20, 0, 1),
    b_lact_mean = rlogis(1, 0, 1),
    b_lact_sd = runif(1, 0, 10),
    b_lact = rlogis(20, 0, 1),
    b_dens_open = rlogis(1, 0, 1),
    b_pdsi_open = rlogis(1, 0, 1),
    b_dens_lact = rlogis(1, 0, 1),
    b_pdsi_lact = rlogis(1, 0, 1)
  )
}

nimble_monitors <- c(
  "b0_mean",
  "b0_sd",
  "b_lact_mean",
  "b_lact_sd",
  "b_dens_open",
  "b_pdsi_open",
  "b_dens_lact",
  "b_pdsi_lact"
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
