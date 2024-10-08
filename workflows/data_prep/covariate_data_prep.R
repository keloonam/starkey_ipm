# Prepares all of the covariates used or potentially used for the original paper

#Packages=======================================================================

require(dplyr); require(readr); require(lubridate); require(tidyr)
require(ggplot2); require(data.table); require(purrr)

#Functions======================================================================

scl <- function(x){
  ((x - mean(x, na.rm = T)) / sd(x, na.rm = T)) %>% return()
}

scl_rb <- function(x, i_range){
  scl_mn <- mean(x[i_range[1]:i_range[2]], na.rm = T)
  scl_sd <-   sd(x[i_range[1]:i_range[2]], na.rm = T)
  out <- (x - scl_mn) / scl_sd
  return(out)
}

#Puma===========================================================================

puma <- tibble(
  yr = 1987:2023,
  pd_reconstruction = NA,
  # pd_mortalities = NA,
  pd_odfw_est = NA,
  pd_logistic = NA,
  pd_full_mean = NA,
  pd_trim_mean = NA
)

pd_reconstruction <- read_csv("data//cougars//cougar_density.csv")

pd_y <- pd_reconstruction$density[1:30]
pd_x <- pd_reconstruction$year[1:30] - 1982
pd_model <- nls(
  pd_y ~ a/(1 + exp(-b * (pd_x - c))), 
  start = list(
    a = 1, 
    b = 0.5, 
    c = 1
  ))

params = coef(pd_model)
pd_logistic <- params[1] / (1 + exp(-params[2] * (1:41 - params[3])))

puma$pd_reconstruction <- pd_reconstruction %>% 
  filter(year > 1986) %>%
  pull(density)
puma$pd_logistic <- pd_logistic[which(pd_reconstruction$year > 1986)]

mortalities <- read_csv("data//cougars//starkey_unit_mort_data.csv") %>%
  mutate(age = case_when(
    Age == 99 ~ NA,
    T ~ Age
  )) %>%
  mutate(dt = mdy(Date)) %>%
  mutate(yr = Year) %>%
  select(yr, dt, age)
puma <- mortalities %>% group_by(yr) %>%
  count() %>%
  mutate(pd_mortalities = n) %>%
  select(-n) %>%
  left_join(puma, .) %>%
  mutate(pd_mortalities = case_when(
    is.na(pd_mortalities) ~ 0,
    T ~ pd_mortalities
  ))

puma$pd_odfw_est[8:37] <- c(
  926,  1059, 1190, 1305, 1396, 1438, 1507, 1593, 1568, 1567, 1585, 1638, 1637, 
  1602, 1589, 1589, 1573, 1534, 1591, 1607, 1681, 1718, 1748, 1745, 1784, 1795, 
  1839, 1891, 1906, 1889
)
# puma$pd_mortalities[1:33] <- c(
#   76,   64,  72,  85,  81,  98,  86,  96,  36,  46,  64,  91, 140, 136, 125, 
#   142, 149, 168, 140, 165, 177, 171, 158, 162, 169, 164, 135,  94, 110, 113, 
#   141, 113, 125
# )

asi <- which(puma$yr %in% 1994:2012)
puma$pd_full_mean <- puma %>%
  mutate(
    recon = scl_rb(pd_reconstruction, c(asi[1], asi[length(asi)])),
    # morts = scl_rb(pd_mortalities, c(asi[1], asi[length(asi)])),
    odfwe = scl_rb(pd_odfw_est, c(asi[1], asi[length(asi)])),
    logis = scl_rb(pd_logistic, c(asi[1], asi[length(asi)]))
  ) %>%
  select(yr, recon, odfwe) %>%
  mutate(pd_full_mean = rowMeans(.[,c(2:3)], na.rm = T)) %>%
  pull(pd_full_mean)
puma$pd_trim_mean <- puma %>%
  mutate(
    recon = scl_rb(pd_reconstruction, c(asi[1], asi[length(asi)])),
    # morts = scl_rb(pd_mortalities, c(asi[1], asi[length(asi)])),
    odfwe = scl_rb(pd_odfw_est, c(asi[1], asi[length(asi)])),
    logis = scl_rb(pd_logistic, c(asi[1], asi[length(asi)]))
  ) %>%
  mutate(recon = case_when(
    yr > 2017 ~ NA,
    T ~ recon
  )) %>%
  select(yr, recon, odfwe) %>%
  mutate(pd_trim_mean = rowMeans(.[,c(2:3)], na.rm = T)) %>%
  pull(pd_trim_mean)
cbp <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
         "#D55E00", "#CC79A7")
puma %>%
  mutate(
    recon = scl_rb(pd_reconstruction, c(asi[1], asi[length(asi)])),
    morts = scl_rb(pd_mortalities, c(asi[1], asi[length(asi)])),
    odfwe = scl_rb(pd_odfw_est, c(asi[1], asi[length(asi)])),
    logis = scl_rb(pd_logistic, c(asi[1], asi[length(asi)])),
    fullm = pd_full_mean,
    trimm = pd_trim_mean
  ) %>%
  select(yr, recon, odfwe, morts, logis, fullm, trimm) %>%
  pivot_longer(cols = 2:7, names_to = "Source", values_to = "Density") %>%
  filter(Source %in% c("fullm", "recon", "odfwe", "logis")) %>%
  filter(!is.na(Density)) %>%
  ggplot(aes(x = yr, y = Density, color = Source, shape = Source)) +
  geom_line() +
  geom_point() +
  xlab("Year") + ylab("N puma (scaled)") + 
  labs(title = "Puma density measures") +
  theme_classic() +
  scale_color_manual(values = cbp) +
  geom_vline(xintercept = 1994, linetype = 2) +
  scale_color_discrete(name = NULL, labels = c("Mean", "Logistic growth", "ODFW Blue Mountain Estimate", "Starkey Reconstruction")) +
  scale_shape_discrete(name = NULL, labels = c("Mean", "Logistic growth", "ODFW Blue Mountain Estimate", "Starkey Reconstruction")) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))
ggsave("figures//puma_covariates.png", dpi = 600, width = 14, height = 9, units = "cm")

#NDVI===========================================================================

avhrr_ndvi <- read_csv("data//climate//avhrr_ndvi_composite.csv") %>%
  mutate(ndvi = (ndvi - 2.261e-5) / 7.548e-5) %>%
  group_by(mn, yr) %>%
  summarise(ndvi_a = mean(ndvi)) %>%
  arrange(yr, mn)
modis_ndvi <- read_csv("data//climate//modis_ndvi_main.csv") %>%
  mutate(ndvi = NDVI) %>%
  group_by(mn, yr) %>%
  summarise(ndvi_m = mean(ndvi)) %>%
  arrange(yr, mn)
ndvi <- modis_ndvi %>% full_join(avhrr_ndvi) %>%
  arrange(yr, mn) %>%
  mutate(ndvi = case_when(
    is.na(ndvi_m) ~ ndvi_a,
    T ~ ndvi_m
  )) %>%
  select(yr, mn, ndvi) %>%
  filter(mn %in% c(5,6,7,8,9)) %>%
  group_by(yr) %>%
  summarise(ndvi = mean(ndvi))

#Prism==========================================================================

prism <- read_csv("data//climate//prism_data_starkey_centroid.csv") %>%
  set_names("dt", "mm_precip", "t_min", "t_mean", "t_max") %>%
  mutate(dt = mdy(dt)) %>%
  mutate(yr = year(dt), mn = month(dt)) %>%
  filter(mn %in% c(5,6,7,8,9)) %>%
  group_by(yr) %>%
  summarise(
    mm_precip = sum(mm_precip),
    temp_mean = mean(t_mean)
  ) %>%
  ungroup() %>%
  # mutate(precip = scl(mm_precip)) %>%
  # mutate(temp = scl(t_mean)) %>%
  select(yr, mm_precip, temp_mean)

# with(prism, plot(t_max ~ t_mean))

#PDI============================================================================

pdi <- read_csv("data//climate//pdi.csv") %>%
  filter(yr >= 1987 & yr < 2024) %>%
  mutate(pdi_growing = (may + june + july + august + september) / 5) %>%
  mutate(pdi_september = september) %>%
  # mutate(pdi_september = scl(september)) %>%
  # mutate(pdi_growing = scl(growing)) %>%
  select(yr, pdi_growing, pdi_september)

#SPEI===========================================================================

spei <- read_csv("data//climate//spei.csv") %>%
  mutate(dt = my(DATA)) %>%
  mutate(yr = year(dt)) %>%
  mutate(mn = month(dt)) %>%
  mutate(spei_1m = SPEI_1) %>%
  mutate(spei_3m = SPEI_3) %>%
  mutate(spei_6m = SPEI_6) %>%
  mutate(spei_12m = SPEI_12) %>%
  mutate(spei_24m = SPEI_24) %>%
  mutate(spei_48m = SPEI_48) %>%
  filter(yr > 1986 & yr < 2024) %>%
  filter(mn == 9) %>%
  select(yr, spei_1m, spei_3m, spei_6m, spei_12m, spei_24m, spei_48m) %>%
  arrange(yr)

#Elk_density====================================================================

est_N <- readRDS("results//submission_2//ipm_1.0_null.Rdata") %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("Ntot", names(.))) %>%
  set_names(as.character(1988:2023)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "est_N") %>%
  group_by(yr) %>%
  summarise(nelk = mean(est_N)) %>%
  mutate(yr = as.numeric(yr))

#Clean-up=======================================================================

covariates <- puma %>%
  full_join(pdi) %>%
  full_join(prism) %>%
  full_join(ndvi) %>%
  full_join(spei) %>%
  full_join(est_N)

saveRDS(covariates, "data//unscaled_covariates.rds")
rm(list = ls())
