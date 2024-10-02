# Build the full IPM data
#Sourcing=======================================================================

source("fb_tests//cjs_data_prep.R")
source("fb_tests//recruitment_data_prep.R")
source("fb_tests//misc_data_prep.R")
climate_cov <- "spei_12m"
puma_cov <- "pd_full_mean"
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
  yr = yr_range,
  fake_data = NA
)

#Load_the_data==================================================================

covariates <- readRDS("data//unscaled_covariates.rds") %>%
  filter(yr %in% yr_range)
misc_data <- readRDS("fb_tests//misc_data.rds")
recruitment <- readRDS("fb_tests//recruitment_data.rds") 
survival <- readRDS("fb_tests//cjs_data.rds")
abundance <- read_csv("data//elk_abundance_ipm_ready.csv") %>%
  filter(year %in% yr_range)
n1_estimates <- readRDS("fb_tests//n1_estimates.rds")
#Pull_and_transform=============================================================

n_added_tibble <- misc_data$n_mov %>%
  filter(year %in% yr_range) %>%
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
n_added_array <- array(0, dim = c(2,2,length(yr_range)))
n_added_array[1,1,] <- n_added_tibble$fe_ca
n_added_array[1,2,] <- n_added_tibble$ma_ca
n_added_array[2,1,] <- n_added_tibble$fe_ad
n_added_array[2,2,] <- n_added_tibble$ma_ad

n_har <- misc_data$n_har %>%
  filter(yr %in% yr_range) %>%
  filter(age == "ad") %>%
  select(yr, sex, harvest) %>%
  pivot_wider(names_from = sex, values_from = harvest) %>%
  full_join(edt) %>% 
  arrange(yr) %>%
  replace_na(list(yr = 0, f = 0, m = 0, fake_data = 0)) %>%
  select(f,m) %>%
  as.matrix() %>% t()

n_min <- array(0, dim = c(2,2,length(yr_range)))
n_min[1,1,] <- misc_data$n_min$calf/2
n_min[1,2,] <- misc_data$n_min$calf/2
n_min[2,1,] <- misc_data$n_min$female
n_min[2,2,] <- misc_data$n_min$male

est_n <- abundance %>%
  mutate(yr = year - (yr_range[1] - 1)) %>%
  mutate(class = case_when(
    age == 1 ~ 1,
    age == 2 & sex == 1 ~ 2,
    age == 2 & sex == 2 ~ 3
  )) %>%
  select(class, yr, mean, tau) %>%
  as.matrix()

cjs_c <- misc_data$count %>%
  mutate(yr = yr - (yr_range[1] - 1)) %>%
  mutate(cnt = count / prob) %>%
  mutate(class = case_when(
    sex == "f" & age == "ad" ~ 2,
    sex == "c" & age == "ca" ~ 1,
    sex == "m" & age == "ad" ~ 3
  )) %>%
  select(class, yr, cnt) %>%
  arrange(yr, class) %>%
  as.matrix()

rec_dt <- matrix(NA, nrow = recruitment$constants$n_years, ncol = 3)
rec_dt[,1] <- recruitment$constants$years - (yr_range[1] - 1)
rec_dt[,2] <- recruitment$data$n_calf
rec_dt[,3] <- recruitment$data$n_cow

est_n1 <- matrix(100, nrow = 2, ncol = 2)
est_n1[1,] <- n1_estimates %>%
  filter(yr == yr_range[1]) %>%
  pull(calf) %>%
  (function(x)x/2)
est_n1[2,1] <- n1_estimates %>%
  filter(yr == yr_range[1]) %>%
  pull(female)
est_n1[2,2] <- n1_estimates %>%
  filter(yr == yr_range[1]) %>%
  pull(male)

clim <- covariates %>%
  filter(yr %in% yr_range) %>%
  pull(climate_cov) %>%
  scl()
puma <- covariates %>%
  filter(yr %in% yr_range) %>%
  pull(puma_cov) %>%
  scl()
nelk <- covariates %>%
  filter(yr %in% yr_range) %>%
  pull(nelk) %>%
  scl()

#Bundle=========================================================================
jags_data <- list(
  nyr = length(yr_range),
  n_added = n_added_array, 
  n_har = n_har,
  n_min = n_min,
  est_n = est_n,
  cjs_c = cjs_c, #[class,year,mean,tau]
  nn = nrow(est_n),
  nc = nrow(cjs_c),
  est_n1 = est_n1,
  nind = survival$constants$nind,
  l = survival$constants$l,
  f = survival$constants$f,
  y = survival$data$y,
  z = survival$data$z,
  m = survival$data$m,
  h = survival$data$h,
  c = survival$data$c,
  rec_dt = rec_dt,
  ny_r = nrow(rec_dt),
  clim = clim,
  puma = puma,
  nelk = nelk
)
