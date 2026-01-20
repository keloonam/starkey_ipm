source("functions//ipm_data_prep_functions.R")

load("data//deprecated//elk_ipm_data_21apr2023.Rdata")
odt <- ipm_data; rm(ipm_data)

# Keep #
ipm_years <- 1988:2023
n_age <- 14

covariates <- readRDS("data//unscaled_covariates.rds")
mscr_rslt <- readRDS("results//mscr_rslt_summ.rds")
recr_rslt <- readRDS("results//recr_rslt_summ.rds")
preg_rslt <- readRDS("results//preg_rslt_summ.rds")
ch_dt <- readRDS("data//capture_handling_data.rds")
mngmnt_dt <- read_csv("data//management_changes.csv") %>%
  filter(year %in% ipm_years) %>%
  mutate(yr = year - 1987)

# Years to filter from each mscr parameter:
fyr_snaf <- c(2017, 2020, 2021, 2023)
fyr_snam <- c(1989, 2020, 2021)
fyr_snca <- c()

#Harvest========================================================================
n_har <- build_age_at_harvest_data(ch_dt, ipm_years, n_age)
haf_i <- apply(n_har[3:n_age,1,], 2, sum) > 0
ham_i <- apply(n_har[3:n_age,2,], 2, sum) > 0
hjf_i <- n_har[2,1,] > 0
hjm_i <- n_har[2,2,] > 0

fyr_shaf <- ipm_years[!haf_i]
fyr_sham <- ipm_years[!ham_i]
fyr_shjf <- ipm_years[!hjf_i]
fyr_shjm <- ipm_years[!hjm_i]

harvest_indicator <- array(1, dim = c(2, 2, dim(n_har)[3]))
harvest_indicator[1,1,] <- harvest_indicator[1,1,] * hjf_i
harvest_indicator[1,2,] <- harvest_indicator[1,2,] * hjm_i
harvest_indicator[2,1,] <- harvest_indicator[2,1,] * haf_i
harvest_indicator[2,2,] <- harvest_indicator[2,2,] * ham_i

#Natural_Survival===============================================================
s_cjs <- mscr_rslt %>%
  ungroup() %>%
  filter(param %in% c("Snaf", "Snam", "Snca")) %>%
  mutate(yr = as.numeric(yr)) %>%
  mutate(year = yr - 1987) %>%
  mutate(age = case_when(
    param == "Snaf" ~ 3,
    param == "Snam" ~ 3,
    param == "Snca" ~ 2
  )) %>%
  mutate(sex = case_when(
    param == "Snaf" ~ 1,
    param == "Snam" ~ 2,
    param == "Snca" ~ 1
  )) %>%
  mutate(rmr = case_when(
    ((param == "Snaf") & (yr %in% fyr_snaf)) ~ T,
    ((param == "Snam") & (yr %in% fyr_snam)) ~ T,
    ((param == "Snca") & (yr %in% fyr_snca)) ~ T,
    T ~ F
  )) %>%
  filter(!rmr) %>%
  select(year, age, sex, mn, tau) %>%
  as.matrix()
#Harvest_Survival===============================================================
h_cjs <- mscr_rslt %>%
  ungroup() %>%
  filter(param %in% c("Shaf", "Sham", "Shjf", "Shjm")) %>%
  mutate(yr = as.numeric(yr)) %>%
  mutate(year = yr - 1987) %>%
  mutate(age = case_when(
    param == "Shaf" ~ 3,
    param == "Sham" ~ 3,
    param == "Shjf" ~ 2,
    param == "Shjm" ~ 2
  )) %>%
  mutate(sex = case_when(
    param == "Shaf" ~ 1,
    param == "Sham" ~ 2,
    param == "Shjf" ~ 1,
    param == "Shjm" ~ 2
  )) %>%
  mutate(rmr = case_when(
    ((param == "Shaf") & (yr %in% fyr_shaf)) ~ T,
    ((param == "Sham") & (yr %in% fyr_sham)) ~ T,
    ((param == "Shjf") & (yr %in% fyr_shjf)) ~ T,
    ((param == "Shjm") & (yr %in% fyr_shjm)) ~ T,
    T ~ F
  )) %>%
  filter(!rmr) %>%
  mutate(mn = mn*-1) %>%
  mutate(age = age - 1) %>%
  select(year, age, sex, mn, tau) %>%
  as.matrix()

#Recruitment====================================================================
r_ratio <- recr_rslt %>%
  filter(yr != 1988) %>%
  filter(yr != 2004) %>%
  mutate(age = NA) %>%
  mutate(sex = NA) %>%
  mutate(year = yr - 1987) %>%
  mutate(tau = 1/sd^2) %>%
  select(year, age, sex, mn, tau) %>%
  as.matrix()
#Pregnancy======================================================================
p_ratio <- preg_rslt %>% as.matrix()
#Management_Movements===========================================================
n_a_mov <- cbind(matrix(0, 2, 1), t(as.matrix(mngmnt_dt[,2:3])))
n_c_mov <- cbind(matrix(0, 2, 1), t(as.matrix(mngmnt_dt[,5:4])))
#Minimums=======================================================================
min_ad <- matrix(0, nrow = 2, ncol = 36)
min_ad[,1:34] <- odt$min_ad
min_ad[1,35:36] <- c(27,22)
min_ad[2,35:36] <- c(9,6)
min_ca <- c(odt$min_ca, 4, 4)

pna <- 0.8^(1:13)/sum(0.8^(1:13))

min_n1 <- matrix(0, nrow = n_age, ncol = 2)
min_n1[1,] <- 166/2
min_n1[2:n_age,1] <- 278 * pna
min_n1[2:n_age,2] <- 32 * pna

est_n1 <- matrix(0, nrow = n_age, ncol = 2)
est_n1[1,] <- 168/2
est_n1[2:n_age,1] <- 404 * pna
est_n1[2:n_age,2] <- 55 * pna

#Puma_Covariates================================================================
pdt <- covariates %>%
  mutate(
    recon = pd_reconstruction,
    orest = pd_odfw_est,
    grwth = pd_logistic
    ) %>%
  select(yr, recon, orest, grwth)
pumas <- tibble(
  yr = pdt$yr,
  reconstruction = (pdt$recon - mean(pdt$recon[8:26])) / sd(pdt$recon[8:26]),
  odfw = (pdt$orest - mean(pdt$orest[8:26])) / sd(pdt$orest[8:26]),
  logistic = (pdt$grwth - mean(pdt$grwth[8:26])) / sd(pdt$grwth[8:26]),
  mn = NA,
  ssd = NA,
  highest = NA,
  lowest = NA,
  mid = NA
) %>%
  filter(yr > 1987)
for(i in 1:nrow(pumas)){
  pumas$mn[i] <- mean(
    c(pumas$reconstruction[i], pumas$odfw[i], pumas$logistic[i]), 
      na.rm = T)
  pumas$ssd[i] <- sd(
    c(pumas$reconstruction[i], pumas$odfw[i], pumas$logistic[i], pumas$mn[i]), 
    na.rm = T)
  pumas$highest[i] <- max(
    c(pumas$reconstruction[i], pumas$odfw[i], pumas$logistic[i], pumas$mn[i]), 
    na.rm = T)
  pumas$lowest[i] <- min(
    c(pumas$reconstruction[i], pumas$odfw[i], pumas$logistic[i], pumas$mn[i]), 
    na.rm = T)
  pumas$mid[i] <- mean(c(pumas$highest[i], pumas$lowest[i]))
}
tmp_ssd_from_mean <- sd(pumas$mn)
tmp_mn <- mean(pumas$mid)
tmp_sd <- sd(pumas$mid)
pumas <- pumas %>%
  mutate(
    recon = scl(reconstruction),
    orest = scl(odfw),
    grwth = scl(logistic),
    mnest = scl(mn),
    mddle = scl(mid),
    sdest = ssd / tmp_ssd_from_mean,
    high = (highest - tmp_mn) / tmp_sd,
    low = (lowest - tmp_mn) / tmp_sd
  ) %>%
  select(yr, recon, orest, grwth, mnest, sdest, high, low)
pumas$orest[1:6] <- min(pumas$orest, na.rm = T)
#Clean-up=======================================================================
fdt <- list(
  s_cjs = s_cjs,
  h_cjs = h_cjs,
  r_ratio = r_ratio,
  p_ratio = p_ratio,
  ns = nrow(s_cjs),
  nh = nrow(h_cjs),
  nr = nrow(r_ratio),
  np = nrow(p_ratio),
  na = n_age,
  n_sight_ca = odt$n_sight_ca,
  n_sight_af = odt$n_sight_af,
  n_sight_am = odt$n_sight_am,
  nn_ca = nrow(odt$n_sight_ca),
  nn_af = nrow(odt$n_sight_af),
  nn_am = nrow(odt$n_sight_am),
  n_a_mov = n_a_mov,
  n_c_mov = n_c_mov,
  n_year = length(1988:2023),
  n_har = n_har,
  min_ad = min_ad,
  min_ca = min_ca,
  af_count = odt$n_f_p_count,
  am_count = odt$n_m_p_count,
  nn_fc = nrow(odt$n_f_p_count),
  nn_mc = nrow(odt$n_m_p_count),
  cg_reconstruction = pumas$recon,
  cg_odfw_estimate = pumas$orest,
  cg_logistic_growth = pumas$grwth,
  cg_mean_estimate = pumas$mnest,
  cg_sd_estiamte = pumas$sdest,
  cg_max_estimate = pumas$high,
  cg_min_estimate = pumas$low,
  pdsi = scl(covariates$pdi_growing[-1]),
  precip = scl(covariates$mm_precip[-1]),
  temp = scl(covariates$temp_mean[-1]),
  spei3 = scl(covariates$spei_3m[-1]),
  spei6 = scl(covariates$spei_6m[-1]),
  spei12 = scl(covariates$spei_12m[-1]),
  ndvi = scl(covariates$ndvi[-1]),
  nelk = c(odt$elk_density, -2.16, -2.16),
  min_n1 = min_n1,
  est_n1 = est_n1,
  harvest_indicator = harvest_indicator
)
saveRDS(fdt, "data//ipm_data_30dec2025.rds")
rm(list = ls())
