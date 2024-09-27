# Build the full IPM data
#Functions======================================================================

sd_to_tau <- function(x){
  return(1/(x^2))
}

scl <- function(x){
  out <- (x - mean(x)) / sd(x)
  return(out)
}

#Empty_tibble===================================================================

edt <- tibble(
  yr = year_range,
  fake_data = NA
)

#Load_the_data==================================================================

covariates <- readRDS("data//unscaled_covariates.rds")
misc_data <- readRDS("data//misc_data.rds")
recruitment <- readRDS("results//recruitment_summary.rds")
survival <- readRDS("results//cjs_summary.rds")
abundance <- read_csv("data//elk_abundance_ipm_ready.csv")

#Pull_and_transform=============================================================

clim <- covariates %>%
  filter(yr %in% year_range) %>%
  pull(climate_cov) %>%
  scl()
puma <- covariates %>%
  filter(yr %in% year_range) %>%
  pull(puma_cov) %>%
  scl()
nelk <- covariates %>%
  filter(yr %in% year_range) %>%
  pull(nelk) %>%
  scl()

n_added_tibble <- misc_data$n_mov %>%
  filter(year %in% year_range) %>%
  mutate(yr = year) %>%
  select(-year) %>%
  full_join(edt) %>% 
  arrange(yr) %>%
  replace_na(list(
    yr = 0, 
    f_ad = 0, 
    m_ad = 0, 
    f_ca = 0, 
    m_ca = 0, 
    fake_data = 0)) %>%
  select(fe_ad, ma_ad, ma_ca, fe_ca)
n_added_array <- array(0, dim = c(2,2,length(year_range)))
n_added_array[1,1,] <- n_added_tibble$fe_ca
n_added_array[1,2,] <- n_added_tibble$ma_ca
n_added_array[2,1,] <- n_added_tibble$fe_ad
n_added_array[2,2,] <- n_added_tibble$ma_ad

n_har <- misc_data$n_har %>%
  filter(yr %in% year_range) %>%
  filter(age == "ad") %>%
  select(yr, sex, harvest) %>%
  pivot_wider(names_from = sex, values_from = harvest) %>%
  full_join(edt) %>% 
  arrange(yr) %>%
  replace_na(list(yr = 0, f = 0, m = 0, fake_data = 0)) %>%
  select(f,m) %>%
  as.matrix() %>% t()

n_min <- array(0, dim = c(2,2,length(year_range)))
n_min[1,1,] <- misc_data$n_min$calf/2
n_min[1,2,] <- misc_data$n_min$calf/2
n_min[2,1,] <- misc_data$n_min$female
n_min[2,2,] <- misc_data$n_min$male

est_n <- abundance %>%
  mutate(yr = year - (year_range[1] - 1)) %>%
  mutate(class = case_when(
    age == 1 ~ 1,
    age == 2 & sex == 1 ~ 2,
    age == 2 & sex == 2 ~ 3
  )) %>%
  select(class, yr, mean, tau) %>%
  as.matrix()

est_bs <- bind_rows(survival) %>%
  filter(yr %in% year_range) %>%
  filter(var %in% c("bsf", "bsc", "bsm")) %>%
  mutate(yr = as.numeric(yr) - (year_range[1] - 1)) %>%
  mutate(class = case_when(
    var == "bsf" ~ 2,
    var == "bsc" ~ 1,
    var == "bsm" ~ 3
  )) %>%
  mutate(tau = sd_to_tau(sd)) %>%
  select(class, yr, mn, tau) %>%
  arrange(yr, class) %>%
  as.matrix()

est_br <- recruitment %>%
  mutate(class = 0) %>%
  mutate(yr = as.numeric(yr) - (year_range[1] - 1)) %>%
  mutate(tau = sd_to_tau(sd)) %>%
  select(class, yr, mn, tau) %>%
  arrange(yr) %>%
  filter(yr > 1) %>%
  as.matrix()

cjs_c <- misc_data$count %>%
  mutate(yr = yr - (year_range[1] - 1)) %>%
  mutate(cnt = count / prob) %>%
  mutate(class = case_when(
    sex == "f" & age == "ad" ~ 2,
    sex == "c" & age == "ca" ~ 1,
    sex == "m" & age == "ad" ~ 3
  )) %>%
  select(class, yr, cnt) %>%
  arrange(yr, class) %>%
  as.matrix()

est_n1 <- matrix(100, nrow = max_age, ncol = 2)
est_n1[1,] <- 168 / 2
est_n1[2:max_age,1] <- 404 / (max_age - 1)
est_n1[2:max_age,2] <- 55 / (max_age - 1)

#Bundle=========================================================================
dtlist <- list(
  nyr = length(year_range),
  max_age = max_age,
  n_added = n_added_array, 
  n_har = n_har,
  n_min = n_min,
  est_n = est_n,
  est_bs = est_bs,
  est_br = est_br,
  cjs_c = cjs_c, #[class,year,mean,tau]
  nn = nrow(est_n),
  ns = nrow(est_bs),
  nr = nrow(est_br),
  nc = nrow(cjs_c),
  clim = clim,
  puma = puma,
  nelk = nelk,
  est_n1 = est_n1
)

#Clean-up=======================================================================
saveRDS(dtlist, "data//ipm_data.rds")
rm(list = ls())
