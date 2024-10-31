require(tidyverse)

load("data//deprecated//elk_ipm_data_21apr2023.Rdata")
odt <- ipm_data; rm(ipm_data)
covariates <- readRDS("data//unscaled_covariates.rds")
misc_data <- readRDS("data//misc_data.rds")
survival <- readRDS("s2//cjs_summary.rds")
abundance <- read_csv("data//elk_abundance_ipm_ready.csv")

scl <- function(x){
  ((x-mean(x))/sd(x)) %>% return()
}

#Survival=======================================================================
s_cjs <- bind_rows(survival[1:3]) %>%
  mutate(year = yr - 1987) %>%
  mutate(age = case_when(
    var == "bsf" ~ 3,
    var == "bsm" ~ 3,
    var == "bsc" ~ 2
  )) %>%
  mutate(sex = case_when(
    var == "bsf" ~ 1,
    var == "bsm" ~ 2,
    var == "bsc" ~ 1
  )) %>%
  select(year, age, sex, mn, tau) %>%
  as.matrix()
#Recruitment====================================================================
r_ratio <- read_rds("s2//recruitment_summary.rds") %>%
  filter(yr != 1988) %>%
  filter(yr != 2004) %>%
  mutate(age = NA) %>%
  mutate(sex = NA) %>%
  mutate(year = yr - 1987) %>%
  select(year, age, sex, mn, tau) %>%
  as.matrix()
#Manegment_Movements============================================================
n_a_mov <- matrix(0, nrow = 2, ncol = 36)
n_c_mov <- matrix(0, nrow = 2, ncol = 36)
n_a_mov[,1:34] <- odt$n_ad_add - abs(odt$n_ad_rem)
n_c_mov[,1:34] <- odt$n_ca_add - abs(odt$n_ca_rem)
n_c_mov[,35] <- -1
#Harvest========================================================================
n_hnt <- array(0, dim = c(4,2,36))
n_hnt[,,1:34] <- odt$n_hnt
n_hnt[4,1,35:36] <- c(4,1)
n_hnt[4,2,35:36] <- c(1,6)
#Minimums=======================================================================
min_ad <- matrix(0, nrow = 2, ncol = 36)
min_ad[,1:34] <- odt$min_ad
min_ad[1,35:36] <- c(27,22)
min_ad[2,35:36] <- c(9,6)
min_ca <- c(odt$min_ca, 4, 4)
#Clean-up=======================================================================
fdt <- list(
  s_cjs = s_cjs,
  r_ratio = r_ratio,
  ns = nrow(s_cjs),
  nr = nrow(r_ratio),
  n_sight_ca = odt$n_sight_ca,
  n_sight_af = odt$n_sight_af,
  n_sight_am = odt$n_sight_am,
  nn_ca = nrow(odt$n_sight_ca),
  nn_af = nrow(odt$n_sight_af),
  nn_am = nrow(odt$n_sight_am),
  n_a_mov = n_a_mov,
  n_c_mov = n_c_mov,
  n_year = length(1988:2023),
  n_har = n_hnt,
  min_ad = min_ad,
  min_ca = min_ca,
  af_count = odt$n_f_p_count,
  am_count = odt$n_m_p_count,
  nn_fc = nrow(odt$n_f_p_count),
  nn_mc = nrow(odt$n_m_p_count),
  cdens = scl(covariates$pd_logistic[-1]),
  pdsi = scl(covariates$pdi_growing[-1]),
  precip = scl(covariates$mm_precip[-1]),
  temp = scl(covariates$temp_mean[-1]),
  spei3 = scl(covariates$spei_3m[-1]),
  spei6 = scl(covariates$spei_6m[-1]),
  spei12 = scl(covariates$spei_12m[-1]),
  ndvi = scl(covariates$ndvi[-1]),
  nelk = c(odt$elk_density, -2.16, -2.16),
  min_n1 = odt$min_n1,
  est_n1 = odt$est_n1
)
saveRDS(fdt, "s2//ipm_data_09oct2024.rds")
rm(list = ls())
