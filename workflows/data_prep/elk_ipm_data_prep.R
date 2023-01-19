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
  "prob_af[33]",
  "prob_am[33]",
  "survival_af[33]",
  "survival_am[33]",
  "survival_ca[33]"
)
cpf_years <- 1:32
cpm_years <- 1:32

#Environment====================================================================

require(tidyverse)

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
  (1 - herd), 2, sum, na.rm = T)[-c(1,34)]
fpc_dat[,4] <- nf.obs / cpf_means

nm.obs <- apply(y * 
  (male) * 
  (1 - (calf == 1)) * 
  (1 - herd), 2, sum, na.rm = T)[-c(1,34)]
mpc_dat[,4] <- nm.obs / cpm_means

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
# stable_density <- mean(cougar_density$density[15:25])
# cougar_density$density[26:nrow(cougar_density)] <- stable_density
# cougar_density_scaled <- as.vector(scale(c(cougar_density$density, stable_density)))

y.t <- cougar_density$density[1:25]
y.ta <- cougar_density$density
x.t <- cougar_density$year[1:25] - 1987
x.ta <- c(cougar_density$year - 1987, 34)
m.cd <- nls(y.t ~ a/(1 + exp(-b * (x.t - c))), start = list(a = 1, 
                                                            b = 0.5, 
                                                            c = 1))

params = coef(m.cd)
y.t2 <- params[1] / (1 + exp(1-params[2] * (x.ta - params[3])))
plot(y.t2, type = "l")
points(y.ta)
cougar_density_scaled <- as.vector(scale(y.t2))

blue_mtns <- tibble(
  year = 1994:2021,
  n = c(926,  1058, 1187, 1301, 1393, 1436, 1502, 1590, 1570, 1570, 1592,
        1646, 1640, 1599, 1592, 1596, 1578, 1541, 1532, 1640, 1703, 1724,
        1748, 1760, 1800, 1807, 1849, 1910),
  Source = "Blue Mountain Estimate"
  ) %>%
  mutate(n = scale(n))

mortalities <- tibble(
  year = 1987:2019,
  n = c(76,   64,  72,  85,  81,  98,  86,  96,  36,  46,  64,  91, 140, 136, 
        125, 142, 149, 168, 140, 165, 177, 171, 158, 162, 169, 164, 135,  94, 
        110, 113, 141, 
        113, 125),
  Source = "Starkey WMU Mortalities"
  ) %>%
  mutate(n = scale(n))

logistic <- tibble(
  year = 1988:2021,
  n = scale(y.t2),
  Source = "Logistic Growth Model"
  )

reconstruction <- tibble(
  year = 1988:2020,
  n = scale(y.ta),
  Source = "Starkey WMU Reconstruction"
)

cougar_tibble <- bind_rows(blue_mtns, mortalities, logistic, reconstruction)
require(ggsci)
ggplot(data = cougar_tibble, aes(x = year, y = n, color = Source, shape = Source, pattern = Source)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = "Year", y = "N Cougars (scaled)", title = "Cougar Abundance") +
  scale_color_jco() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = T))
ggsave("cougar_density_plot.png", width = 5, height = 3, units = "in", dpi = 300)

# plot(cougar_density_scaled, type = "l", xlim = c(0, 35), ylim = c(-2,1.5))
# points(scale(y.ta), col = "blue")
# points(scale(c(rep(NA, 7), n_cougars)), col = "red")

#Density========================================================================

load("results//ipm_result_11oct2022_R_null.Rdata")
scaled_density <- as.numeric(scale(summary(rslt)$statistics[103:136,1]))
scaled_N_AF <- as.numeric(scale(summary(rslt)$statistics[35:68,1]))
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

#Combine========================================================================

ipm_data <- list(
  s_cjs = cjs_dat,
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
  summer_temp = sum_temp,
  winter_temp = win_temp,
  summer_precip = summer_precip,
  august_precip = august_precip,
  winter_precip = winter_precip,
  n_f_p_count = fpc_dat,
  n_m_p_count = mpc_dat,
  ratio_counts = ratio_counts,
  cougar_density = cougar_density_scaled,
  elk_density = scaled_density,
  af_density = scaled_N_AF,
  palmer_index = pdi,
  ndvi_avhrr = ndvi_avhrr,
  ndvi_modis = ndvi_modis
)

save(ipm_data, file = "data//elk_ipm_data_05jan2023.Rdata")

