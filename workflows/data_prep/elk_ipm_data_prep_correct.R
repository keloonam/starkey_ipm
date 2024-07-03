# IPM data prep -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# User Specified
n_yrs <- 34
cjs_yrs <- c(min_yr = 1, max_yr = 33)
ratio_yrs <- c(min_yr = 2, max_yr = 34)
first_year <- 1988
start_year <- first_year
end_year <- 2021
yr_a <- first_year
yr_z <- end_year

cjs_file <- "results//survival//cjs_rslt_13sep2022.Rdata"
ratio_file <- "results//recruitment_years_2to31.Rdata"
abundance_file <- "data//elk_abundance_ipm_ready.csv"

# Survival estimates to remove from CJS
# Parameters either did not converge or 
# are from years with unrecorded id's removed
cjs_removals <- c(
  "prob_af[34]",
  "prob_am[34]",
  "survival_af[34]",
  "survival_am[34]",
  "survival_ca[34]",
  "survival_af[32]",
  "survival_am[32]",
  "survival_ca[32]",
  "survival_af[29]"
)
cpf_years <- 1:33
cpm_years <- 1:33

#Environment====================================================================

require(tidyverse)

scale_clean <- function(x){
  out <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  return(out)
}

logit <- function(x){
  out <- log(x / (1 - x))
  return(out)
}

calculate_tau <- function(x){
  out <- 1/sd(x)^2
  return(out)
}

#CJS_Survival===================================================================

load(cjs_file)
cjs_raw <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(starts_with("survival_")) %>%
  select(-contains(cjs_removals))
cp_raw <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(starts_with("prob_")) %>%
  select(-contains(cjs_removals))
rm(rslt)

# first year is survival to 1989
cjs_years <- as.numeric(str_extract_all(names(cjs_raw), "[0-9]+")) + 1
cjs_males <- grep("am", names(cjs_raw))
cjs_calfs <- grep("ca", names(cjs_raw))
cjs_means <- unlist(map(cjs_raw, mean))
cjs_sd    <- unlist(map(cjs_raw, sd))

cjs_dat <- matrix(1, nrow = length(cjs_years), ncol = 5)
cjs_dat[,1] <- cjs_years
cjs_dat[,2] <- 3
cjs_dat[cjs_calfs,2] <- 2
cjs_dat[cjs_males,3] <- 2
cjs_dat[,4] <- cjs_means
cjs_dat[,5] <- 1/(cjs_sd^2)

cp_years <- as.numeric(unlist(str_extract_all(names(cp_raw), "[0-9]+"))) + 1
cp_years_m <- cp_years[max(cpf_years) + cpm_years]
cp_years_f <- cp_years[cpf_years]

cp_means <- unlist(map(cp_raw, mean))
cpf_means <- cp_means[cpf_years]
cpm_means <- cp_means[max(cpf_years) + cpm_years]

fpc_dat <- matrix(nrow = length(cpf_means), ncol = 5)
mpc_dat <- matrix(nrow = length(cpm_means), ncol = 5)
fpc_dat[,1] <- cp_years_f
mpc_dat[,1] <- cp_years_m
p_cjs <- matrix(NA, nrow = length(cp_years), ncol = 8)
p_cjs[,1] <- cp_years 
p_cjs[,3] <- c(rep(1, length(cp_years_f)), rep(2, length(cp_years_m)))
p_cjs[,4] <- cp_raw %>% 
  map(logit) %>%
  map(mean) %>%
  unlist()
p_cjs[,5] <- cp_raw %>% 
  map(logit) %>%
  map(calculate_tau) %>%
  unlist()
p_cjs[,7] <- cp_raw %>%
  map(logit) %>%
  map(quantile, 0.025) %>%
  unlist()
p_cjs[,8] <- cp_raw %>% 
  map(logit) %>%
  map(quantile, 0.975) %>%
  unlist()

load("data//elk_data.Rdata")
y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
male <- elk_data$sex_tib %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)
calf <- elk_data$age_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()
rm(elk_data)

nf.obs <- apply(y * 
  (1 - male) * 
  (1 - (calf == 1)) * 
  (1 - herd), 2, sum, na.rm = T)[-c(1)]
fpc_dat[,4] <- nf.obs / cpf_means

nm.obs <- apply(y * 
  (male) * 
  (1 - (calf == 1)) * 
  (1 - herd), 2, sum, na.rm = T)[-c(1)]
mpc_dat[,4] <- nm.obs / cpm_means

p_cjs[,6] <- c(nf.obs, nm.obs)

#CJS data prep==================================================================

load("data//elk_data.Rdata")

y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  as.matrix()

l <- elk_data$hnt_tib %>%
  arrange(id) %>%
  mutate('2021' = 1) %>%
  select(as.character(yr_a:yr_z)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  as.matrix() 

y <- w + y # harvests count as observed alive that year (did not die naturally)

male_id <- elk_data$sex_tib %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)

calf <- elk_data$age_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  as.matrix()

herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(yr_a:yr_z)) %>%
  mutate(across(all_of(as.character(yr_a:yr_z)), ~ (.x != "main")*1)) %>%
  as.matrix()

gone_elk <- apply(herd, 1, sum, na.rm = T) == ncol(herd)
unseen_elk <- apply(y, 1, sum) == 0
f <- apply(y, 1, function(x) min(which(x != 0)))
weird_elk <- f == l
rm_elk <- which((gone_elk + unseen_elk + weird_elk) != 0)

y <- y[-rm_elk,]
l <- l[-rm_elk]
w <- w[-rm_elk,]
male_id <- male_id[-rm_elk]
calf <- calf[-rm_elk,]
herd <- herd[-rm_elk,]
f <- f[-rm_elk]
calf <- (calf == 1) * 1
calf[is.na(calf)] <- 0

female <- matrix(0, nrow = nrow(calf), ncol = ncol(calf))
male   <- matrix(0, nrow = nrow(calf), ncol = ncol(calf))
for(i in 1:length(male_id)){
  if(male_id[i] == 1){
    male[i,] <- (1 - calf[i,])
  }else{
    female[i,] <- (1 - calf[i,])
  }
}

z <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
l_k <- rep(NA, nrow(y))
for(i in 1:nrow(y)){
  l_k[i] <- max(which(y[i,] == 1))
  z[i, (f[i]):l_k[i]] <- 1
}

#Ratio_Recruitment==============================================================

load(ratio_file)
ratio_raw <- rslt$BUGSoutput$summary
rm(rslt)

# years that recruitment model didn't converge (no data)
ratio_removals <- c(3,4,5,6,7,9,16,33)

ratio_data <- ratio_raw[grep("R", dimnames(ratio_raw)[[1]]),]

r_ratio <- matrix(NA, nrow = nrow(ratio_data), ncol = 5)
r_ratio[,1] <- ratio_yrs[1]:ratio_yrs[2]
r_ratio[,4] <- ratio_data[,1]
r_ratio[,5] <- 1/ratio_data[,2]^2
r_ratio <- r_ratio[-ratio_removals,]

r_data <- read_csv("data//min_n_handle_summaries.csv") %>%
  mutate(yr = year - yr_a + 1) %>%
  select(yr, ca, fe_ad) %>%
  as.matrix() 
r_data[33,2] <- 13 # One year there were more calves than cows - fix for binom

#Abundance======================================================================

n_raw <- read_csv(abundance_file)

n_sight_ca <- n_raw %>%
  filter(age == 1) %>%
  as.matrix(.)
n_sight_ca[,1] <- n_sight_ca[,1] - first_year + 1

n_sight_af <- n_raw %>%
  filter(age == 2 & sex == 1) %>%
  as.matrix(.)
n_sight_af[,1] <- n_sight_af[,1] - first_year + 1

n_sight_am <- n_raw %>%
  filter(age == 2 & sex == 2) %>%
  as.matrix(.)
n_sight_am[,1] <- n_sight_am[,1] - first_year + 1

#Min_N==========================================================================

min_dat <- read_csv("data//min_n_handle_summaries.csv") %>%
  pivot_longer(cols = c(fe_ad, ma_ad, ca)) %>%
  mutate(year = year - 1987)

n_fg_c <- matrix(NA, nrow = sum(min_dat$name == "ca"), ncol = 4)
n_fg_c[,1] <- min_dat$year[min_dat$name == "ca"]
n_fg_c[,4] <- min_dat$value[min_dat$name == "ca"]

n_fg_m <- matrix(NA, nrow = sum(min_dat$name == "ma_ad"), ncol = 4)
n_fg_m[,1] <- min_dat$year[min_dat$name == "ma_ad"]
n_fg_m[,4] <- min_dat$value[min_dat$name == "ma_ad"]

n_fg_f <- matrix(NA, nrow = sum(min_dat$name == "fe_ad"), ncol = 4)
n_fg_f[,1] <- min_dat$year[min_dat$name == "fe_ad"]
n_fg_f[,4] <- min_dat$value[min_dat$name == "fe_ad"]

ratio_counts <- matrix(NA, nrow = nrow(n_fg_c), ncol = 3)
ratio_counts[,1] <- n_fg_c[,1]
ratio_counts[,2] <- n_fg_f[,4]
ratio_counts[,3] <- n_fg_c[,4]

#Harvested======================================================================

load("data//elk_harvest_data.Rdata")
n_hnt <- array(data = 0, dim = c(4,2,n_yrs))
hntdat <- as.matrix(hntdat)
for(i in 1:(nrow(hntdat))){
  n_hnt[hntdat[i,2], hntdat[i,3], hntdat[i,1]] <- hntdat[i,4]
}

#Animal_movement================================================================

# Build 3 dim n_mov array [age,sex,year]
mov_dat <- read_csv("data//mov_data_handle_summaries.csv") %>%
  as.matrix()
mov_dat[,1] <- mov_dat[,1] - 1987
n_mov <- array(0, dim = c(3,2,n_yrs))
for(i in 1:nrow(mov_dat)){
  n_mov[1,1,mov_dat[i,1]] <- mov_dat[i,5]
  n_mov[1,2,mov_dat[i,1]] <- mov_dat[i,4]
  n_mov[3,1,mov_dat[i,1]] <- mov_dat[i,2]
  n_mov[3,2,mov_dat[i,1]] <- mov_dat[i,3]
}

n_add <- array(0, dim = dim(n_mov))
n_rem <- array(0, dim = dim(n_mov))

for(i in 1:dim(n_mov)[1]){
  for(j in 1:dim(n_mov)[2]){
    for(k in 1:dim(n_mov)[3]){
      if(n_mov[i,j,k] > 0){
        n_add[i,j,k] <- n_mov[i,j,k]
      }
      if(n_mov[i,j,k] < 0){
        n_rem[i,j,k] <- n_mov[i,j,k]
      }
    }
  }
}

n_ad_rem <- abs(n_rem[2,,] + n_rem[3,,])
n_ca_rem <- abs(n_rem[1,,])
n_ad_add <- n_add[2,,] + n_add[3,,]
n_ca_add <- n_add[1,,]

#Min_N_alive_capture_database===================================================

load("data//elk_minimum_count_data.Rdata")

min_counts$af[1] <- 278
min_counts$am[1] <- 32
min_counts$ca[1] <- 166

min_ad <- matrix(0, nrow = 2, ncol = n_yrs)
min_ca <- rep(0, n_yrs)
for(i in 1:nrow(min_counts)){
  j <- min_counts$year[i] - 1987
  min_ca[j] <- min_counts$ca[i]
  min_ad[1,j] <- min_counts$af[i]
  min_ad[2,j] <- min_counts$am[i]
}

#Climate========================================================================

monthly_temp <- read_csv("data//climate//ElkClimateChange_microclimate.csv") %>%
  mutate(add_year = case_when(
    Month >= 11 ~ Year + 1,
    T ~ Year
  ))

mean_ann_temp <- monthly_temp %>%
  group_by(Year) %>%
  summarise(temp = mean(Starkey_HQ_temp_C)) %>%
  mutate(temp = scale(temp))
ann_temp <- rep(0, n_yrs)
for(i in 1:nrow(mean_ann_temp)){
  ann_temp[mean_ann_temp$Year[i] - 1987] <- mean_ann_temp$temp[i]
}

summer_temp <- monthly_temp %>%
  filter(between(Month, 7, 9)) %>%
  group_by(Year) %>%
  summarise(st = mean(Starkey_HQ_temp_C)) %>%
  mutate(st = scale(st))
sum_temp <- rep(0, n_yrs)
for(i in 1:nrow(summer_temp)){
  sum_temp[summer_temp$Year[i] - 1987] <- summer_temp$st[i]
}

winter_temp <- monthly_temp %>%
  filter(!between(Month, 3, 10)) %>%
  group_by(Year) %>%
  summarise(wt = mean(Starkey_HQ_temp_C)) %>%
  mutate(wt = scale(wt))
win_temp <- rep(0, n_yrs)
for(i in 1:nrow(winter_temp)){
  win_temp[winter_temp$Year[i] - 1987] <- winter_temp$wt[i]
}



get_precip <- function(climate_df, season_months, n_yrs){
  seasonal_precip <- climate_df %>%
    filter(Month %in% season_months) %>%
    group_by(add_year) %>%
    summarise(seasonal_precip = sum(SNOTEL_CL_precipt_mm)) %>%
    mutate(seasonal_precip = scale(seasonal_precip))
  out <- rep(0, n_yrs)
  for(i in 1:nrow(seasonal_precip)){
    out[seasonal_precip$add_year[i] - 1987] <- seasonal_precip$seasonal_precip[i]
  }
  return(out)
}

summer_precip <- get_precip(
  climate_df = monthly_temp, 
  season_months = c(7, 8, 9),
  n_yrs = n_yrs
)
august_precip <- get_precip(
  climate_df = monthly_temp, 
  season_months = c(8),
  n_yrs = n_yrs
)
winter_precip <- get_precip(
  climate_df = monthly_temp, 
  season_months = c(1, 2, 12),
  n_yrs = n_yrs
)

#Cougars========================================================================

cougar_density <- read_csv("data//cougars//cougar_density.csv")

y.t <- cougar_density$density[1:25]
y.ta <- cougar_density$density
x.t <- cougar_density$year[1:25] - 1987
x.ta <- c(cougar_density$year - 1987, 34)
m.cd <- nls(y.t ~ a/(1 + exp(-b * (x.t - c))), start = list(a = 1, 
                                                            b = 0.5, 
                                                            c = 1))

params = coef(m.cd)
y.t2 <- params[1] / (1 + exp(1-params[2] * (x.ta - params[3])))
puma_derived <- as.vector(scale(y.t2))

reconstruction <- read_csv("data//cougars//cougar_density.csv") %>%
  filter(year < 2016) %>%
  mutate(n = scale_clean(density)) %>%
  mutate(method = "reconstruction") %>%
  select(year, method, n)

blue_mtns <- tibble(
  year = 1994:2021,
  n = c(926,  1058, 1187, 1301, 1393, 1436, 1502, 1590, 1570, 1570, 1592,
        1646, 1640, 1599, 1592, 1596, 1578, 1541, 1532, 1640, 1703, 1724,
        1748, 1760, 1800, 1807, 1849, 1910),
  method = "odfw_estimate"
) %>%
  mutate(n = scale_clean(n))

mortalities <- tibble(
  year = 1987:2019,
  n = c(76,   64,  72,  85,  81,  98,  86,  96,  36,  46,  64,  91, 140, 136, 
        125, 142, 149, 168, 140, 165, 177, 171, 158, 162, 169, 164, 135,  94, 
        110, 113, 141, 
        113, 125),
  method = "mortalities"
) %>%
  mutate(n = scale_clean(n))

puma_data <- bind_rows(reconstruction, blue_mtns, mortalities)
full_puma_data <- puma_data %>%
  group_by(year) %>%
  summarise(n = mean(n)) %>%
  ungroup() %>%
  mutate(method = "group_mean") %>%
  bind_rows(puma_data) %>%
  add_row(year = 1988:2021, method = "logistic_regression", n = puma_derived)

puma_mean <- full_puma_data %>% 
  filter(method == "group_mean") %>%
  pull(n) %>%
  scale_clean()

#Density========================================================================

load("results//ipm//ipm_result_14mar2023_R_null.Rdata")
scaled_density <- as.numeric(scale(summary(rslt)$statistics[104:137,1]))
scaled_N_AF <- as.numeric(scale(summary(rslt)$statistics[36:69,1]))
rm(rslt)

#Palmer_Drought_Index===========================================================

pdi_full <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  filter(year >= 1988) %>%
  pivot_longer(2:13, names_to = "month")

pdi <- pdi_full %>%
  filter(month %in% c("september")) %>%
  select(value) %>%
  scale() %>%
  as.vector()

pdi <- pdi[-35]

#NDVI===========================================================================

avhrr <- read_csv("data//climate//ndvi_avhrr.csv") %>%
  mutate(source = "avhrr") %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI)) %>%
  filter(yr < 2015)

modis <- read_csv("data//climate//ndvi_modis.csv") %>%
  mutate(source = "modis") %>%
  filter(SummaryQA <= 1) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(NDVI))

combo <- inner_join(avhrr, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  arrange(yr, mn) %>%
  filter(yr < 2015) %>%
  filter(mn > 4, mn < 10)

m1 <- with(combo, lm(ndvi_m ~ ndvi_a))
avhrr_m1 <- avhrr %>%
  mutate(
    ndvi = ndvi*m1$coefficients[2] + m1$coefficients[1]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))
ndvi_modis <- full_join(avhrr_m1, modis, by = c("mn", "yr"), suffix = c("_a", "_m")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_m) ~ ndvi_m,
    T ~ ndvi_a
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi)) %>%
  filter(yr %in% 1988:2021) %>%
  select(summer_ndvi) %>%
  scale() %>%
  as.vector()

m2 <- with(combo, lm(ndvi_a ~ ndvi_m))
modis_m2 <- modis %>%
  mutate(
    ndvi = ndvi*m2$coefficients[2] + m2$coefficients[1]
  ) %>%
  group_by(mn, yr) %>%
  summarise(ndvi = mean(ndvi))
ndvi_avhrr <- full_join(modis_m2, avhrr, by = c("mn", "yr"), suffix = c("_m", "_a")) %>%
  mutate(ndvi = case_when(
    !is.na(ndvi_a) ~ ndvi_a,
    T ~ ndvi_m
  )) %>%
  filter(mn > 4, mn < 10) %>%
  group_by(yr) %>%
  summarise(summer_ndvi = mean(ndvi)) %>%
  filter(yr %in% 1988:2021) %>%
  select(summer_ndvi) %>%
  scale() %>%
  as.vector()

#PRISM_data=====================================================================

prism <- read_csv("data/climate/prism_data_starkey.csv", skip = 10) %>%
  data.table::setnames(
    old = names(.),  
    new = c("date", "ppt_inches", "tmin_f", "tmean_f", "tmax_f", "vpdmin", 
            "vpdmax")) %>%
  mutate(precip = ppt_inches * 25.4) %>%
  mutate(temp = (tmean_f - 32) / 1.8) %>%
  separate(date, into = c("year", "month"), sep = "-") %>%
  mutate(year = as.numeric(year)) %>%
  mutate(month = as.numeric(month)) %>%
  select(year, month, precip, temp) %>%
  filter(year > 1987 & year < 2022) %>%
  filter(month %in% c(7, 8, 9)) %>%
  group_by(year) %>%
  summarise(summer_precip = sum(precip), summer_temp = mean(temp))


#Adjust min N===================================================================

min_n1 <- matrix(0, nrow = 4, ncol = 2)
est_n1 <- matrix(0, nrow = 4, ncol = 2)

min_n1[1,1] <- min_ca[1] / 2
min_n1[1,2] <- min_ca[1] / 2

min_n1[2,1] <- min_ad[1,1] / 3
min_n1[3,1] <- min_ad[1,1] / 3
min_n1[4,1] <- min_ad[1,1] / 3

min_n1[2,2] <- min_ad[2,1] / 3
min_n1[3,2] <- min_ad[2,1] / 3
min_n1[4,2] <- min_ad[2,1] / 3

est_n1[1,1] <- 168 / 2
est_n1[1,2] <- 168 / 2

est_n1[2,1] <- 404 / 3
est_n1[3,1] <- 404 / 3
est_n1[4,1] <- 404 / 3

est_n1[2,2] <- 55 / 3
est_n1[3,2] <- 55 / 3
est_n1[4,2] <- 55 / 3

#Combine========================================================================



ipm_data <- list(
  y = y,
  z = z,
  male = male,
  female = female,
  calf = calf,
  herd = herd,
  l = l,
  f = f,
  n_ind = nrow(y),
  cjs_n_male = nf.obs,
  cjs_n_female = nm.obs,
  s_cjs = cjs_dat,
  p_cjs = p_cjs,
  r_ratio = r_ratio,
  n_sight_ca = n_sight_ca,
  n_sight_af = n_sight_af,
  n_sight_am = n_sight_am,
  n_fg_c = n_fg_c,
  n_fg_m = n_fg_m,
  n_fg_f = n_fg_f,
  n_ad_add = n_ad_add,
  n_ad_rem = n_ad_rem,
  n_ca_add = n_ca_add,
  n_ca_rem = n_ca_rem,
  n_hnt = n_hnt,
  min_ca = min_ca,
  min_ad = min_ad,
  annual_temp = ann_temp,
  summer_temp = prism$summer_temp,
  winter_temp = win_temp,
  summer_precip = prism$summer_precip,
  august_precip = august_precip,
  winter_precip = winter_precip,
  n_f_p_count = fpc_dat,
  n_m_p_count = mpc_dat,
  ratio_counts = ratio_counts,
  puma_derived = puma_derived,
  puma_mean = puma_mean[-length(puma_mean)],
  puma_reconstruction = c(scale_clean(cougar_density$density), 0),
  elk_density = scaled_density,
  af_density = scaled_N_AF,
  palmer_index = pdi,
  ndvi_avhrr = ndvi_avhrr,
  ndvi_modis = ndvi_modis,
  min_n1 = min_n1,
  est_n1 = est_n1,
  r_data = r_data[-33,] # this year has more calves than cows
)

save(ipm_data, file = "data//elk_ipm_data_14jun2024.Rdata")
rm(list = ls())
