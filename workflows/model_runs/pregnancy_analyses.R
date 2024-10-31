# Pregnancy analysis
# Kenneth Loonam
# May 2023

#Variables======================================================================


niter <- 50000
nburn <- 10000
nthin <- 100

chain <- 3
save_file <- "results//pregnancy_analysis_13oct2024.rds"

#Packages=======================================================================

require(tidyverse); require(nimble)

#Load stuff=====================================================================
rd <- read.csv("data//pregnancy_rates.csv") %>%
  as_tibble()

ipm_data <- readRDS("s2//ipm_data_25sep2024.rds")
scl <- function(x){
  (x-mean(x)) / sd(x)
}
source("models//misc//pregnancy.R")
#Prep data======================================================================
cov_dat <- tibble(
  pdi = ipm_data$spei12[2:34],
  elk = ipm_data$nelk[2:34],
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
  "a_preg_lac_prime",
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
  "bpdi_lac_old",
  "p",
  "bpdi_lac_old",
  # "rsdl",
  "preg_young"
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
  nchains = chain,
  thin = nthin
)
  
mcmcplots::mcmcplot(nimble_results)

saveRDS(nimble_results, file = save_file)
nimble_results <- readRDS(save_file)
#Clean results==================================================================

row_names <- dimnames(nimble_results[[1]])[[2]]
row_names %>% grepl("p\\[", .) %>% sum()

row_keep <- rd %>%
  mutate(keep = (lactating == T) * (age_class == "prime")) %>%
  pull(keep)

preg_years_to_pull <- rd %>% 
  filter(row_keep == T) %>%
  left_join(., cov_dat) %>% pull(yr)
nelk_from_preg_years <- n_af %>%
  filter(year %in% preg_years_to_pull) %>%
  pull(mean_density)
point_data <- rd %>% 
  mutate(yr = yr - 1988) %>%
  filter(row_keep == T) %>%
  left_join(., cov_dat) %>%
  mutate(pdi = pdi * 1.054813 - 0.5071725) %>%
  mutate(elk = nelk_from_preg_years)

p_point_dat <- nimble_results %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("a_preg_lac_prime", names(.))) %>%
  map(quantile, probs = c(.025, .5, .975)) %>%
  bind_rows() %>%
  set_names("lci", "mci", "uci") %>%
  mutate(yr = 1989:2021) %>%
  mutate(spei = cov_dat$pdi) %>%
  mutate(nelk = cov_dat$elk)

ggplot(p_point_dat, aes(x = nelk, y = mci)) +
  geom_pointrange(aes(ymax = uci, ymin = lci))

#Clean results==================================================================

row_names <- dimnames(nimble_results[[1]])[[2]]

cov_dat_problem_child <- nimble_results %>%
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

# ggplot(cov_dat, aes(x = covariate, y = med, color = age, shape = lactating)) +
#   geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = .5)) +
#   theme(legend.title = element_blank()) +
#   labs(x = "", y = "Covariate value", title = "Effects on pregnancy rates") +
#   theme_classic() +
#   geom_hline(yintercept = 0, linetype = "dotted")
#   
# source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")
# 
# ggplot(prg_dat, aes(x = `Age/lactation class`, y = val, color = `Age/lactation class`, fill = `Age/lactation class`)) +
#   geom_flat_violin() +
#   theme_classic() +
#   theme(legend.position = "none") +
#   labs(x = "", y = "Posterior distribution", title = "Mean pregnancy rates") +
#   coord_flip()
# 
# samp_size <- rd %>%
#   mutate(ratio = paste0(
#     " ",
#     round(pregnancy_rate * n_observations),
#     "/",
#     n_observations
#     )) %>%
#   mutate(class = case_when(
#     age_class == "old" & lactating == TRUE ~ "Old - lactating",
#     age_class == "old" & lactating == FALSE ~ "Old",
#     age_class == "prime" & lactating == TRUE ~ "Prime - lactating",
#     age_class == "prime" & lactating == FALSE ~ "Prime",
#     age_class == "young" & lactating == FALSE ~ "Young",
#     age_class == "young" & lactating == TRUE ~ "Young - lactating"
#   )) %>%
#   select(yr, ratio, class) %>%
#   pivot_wider(names_from = class, values_from = ratio) %>%
#   write.csv("data//pregnancy_sample_sizes.csv")


preg <- readRDS(save_file)

p_eff_post <- nimble_results %>%
  map(as_tibble) %>%
  bind_rows() 

# %>%
#   pivot_longer(1:ncol(.), names_to = "parameter") %>%
#   filter(parameter %in% c(
#     "b0_mean", 
#     "b_lac_prm_mean", "bden_lac_prm", "bpdi_lac_prm",
#     "b_dry_prm_mean", "bden_dry_prm", "bpdi_dry_prm",
#     "b_lac_old_mean", "bden_lac_old", "bpdi_lac_old",
#     "b_dry_old_mean", "bden_dry_old", "bpdi_dry_old",
#     "b_dry_yng_mean", "bden_dry_yng", "bpdi_dry_yng"
#   )) %>%
#   mutate(id = sort(rep(1:(nrow(.)/4), 4))) %>%
#   group_by(parameter) %>%
#   summarise(
#     lci = quantile(value, .025),
#     mci = quantile(value, .5),
#     uci = quantile(value, .975),
#     pod = mean(value > 0)
#   )


expit <- function(x){
  1/(1+exp(-x))
}

#####
x <- seq(-3, 3, length.out = 1000)

# PDSI t marginal plot data
p_pt_line <- matrix(data = NA, nrow = length(x), ncol = 4)
p_pt_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(p_eff_post$b0_mean + x[i]*p_eff_post$bpdi_lac_prm + p_eff_post$b_lac_prm_mean) 
  p_pt_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
p_pt_line <- as_tibble(p_pt_line)
names(p_pt_line) <- c("val", "lci", "mci", "uci")
p_pt_line$val <- p_pt_line$val * sd(pdi_dat$pdi) + mean(pdi_dat$pdi)

# PDSI t marginal plot data
p_ed_line <- matrix(data = NA, nrow = length(x), ncol = 4)
p_ed_line[,1] <- x
for(i in 1:length(x)){
  y <- expit(p_eff_post$b0_mean + x[i]*p_eff_post$bden_lac_prm + p_eff_post$b_lac_prm_mean) 
  p_ed_line[i,2:4] <- quantile(y, probs = c(0.025, 0.5, 0.975))
}
p_ed_line <- as_tibble(p_ed_line)
names(p_ed_line) <- c("val", "lci", "mci", "uci")
p_ed_line$val <- p_ed_line$val * sd(n_af$mean_density) + mean(n_af$mean_density)

ggplot(cov_dat_problem_child, aes(x = covariate, y = med, color = age, shape = lactating)) +
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

#Writing and need results summary===============================================

prg_dat %>%
  group_by(`Age/lactation class`) %>%
  summarise(
    lci = quantile(val, .025),
    mci = quantile(val, 0.50),
    uci = quantile(val, .975)
    )
np_cov_dat <- preg_rd %>%
  filter(lactating == TRUE) %>%
  filter(age_class == "prime") %>%
  full_join(n_af) %>% 
  select(yr, lactating, pregnancy_rate, n_observations, age_class, mean_density) %>%
  filter(!is.na(age_class)) %>%
  full_join(p_dat)
