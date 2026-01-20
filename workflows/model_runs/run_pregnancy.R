require(tidyverse); require(nimble)
rd <- read.csv("data//pregnancy_rates.csv") %>%
  as_tibble() %>%
  mutate(preg = round(pregnancy_rate * n_observations)) %>%
  group_by(yr, age_class) %>%
  summarise(no = sum(n_observations), np = sum(preg))

nimble_data <- list(
  nobs = rd$no,
  preg = rd$np
)  

nimble_constants <- list(
  nr = nrow(rd)
)

nimble_inits <- function(){
  list(
    p = runif(nrow(rd), 0, 1)
  )
}

nimble_monitors <- c("p")
nimble_code <- nimble::nimbleCode({
  # GLM
  for(i in 1:nr){
    # Prior
    p[i] ~ dunif(0, 1)
    # Model  
    preg[i] ~ dbin(p[i], nobs[i])
  }
}
)
#Run model======================================================================

nimble_results <- nimbleMCMC(
  code = nimble_code,
  constants = nimble_constants,
  data = nimble_data,
  inits = nimble_inits,
  monitors = nimble_monitors,
  niter = 25000,
  nburnin = 5000,
  nchains = 3,
  thin = 10,
  progressBar = T,
  summary     = F,
  check       = F
)

# mcmcplots::mcmcplot(nimble_results)

saveRDS(nimble_results, file = "results//preg_rslt.rds")

nimble_results %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pivot_longer(cols = 1:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(
    yr = rep(rd$yr, 6000),
    age_ch = rep(rd$age_class, 6000)
    ) %>%
  mutate(bv = logit(val)) %>%
  mutate(year = yr - 1986) %>%
  mutate(age = case_when(
    age_ch == "old" ~ 3,
    age_ch == "prime" ~ 2,
    age_ch == "young" ~ 1
  )) %>%
  mutate(sex = NA) %>% 
  group_by(year, age, sex) %>%
  summarise(mn = mean(bv), tau = 1/sd(bv)^2) %>%
  saveRDS(file = "results//preg_rslt_summ.rds")

