# Prepares all of the covariates used or potentially used for the original paper

#Packages=======================================================================

require(dplyr); require(readr); require(lubridate); require(tidyr)
require(ggplot2); require(data.table); require(purrr)

#Functions======================================================================

scl <- function(x){
  ((x - mean(x)) / sd(x)) %>% return()
}

scl_rb <- function(x, i_range){
  scl_mn <- mean(x[i_range[1]:i_range[2]])
  scl_sd <-   sd(x[i_range[1]:i_range[2]])
  out <- (x - scl_mn) / scl_sd
  return(out)
}

#Puma===========================================================================

puma <- tibble(
  yr = 1987:2023,
  pd_reconstruction = NA,
  pd_mortalities = NA,
  pd_odfw_est = NA,
  pd_logistic = NA
)

pd_reconstruction <- read_csv("data//cougars//cougar_density.csv")

pd_y <- pd_reconstruction$density[1:25]
pd_x <- pd_reconstruction$year[1:25] - 1987
pd_model <- nls(
  pd_y ~ a/(1 + exp(-b * (pd_x - c))), 
  start = list(
    a = 1, 
    b = 0.5, 
    c = 1
  ))

params = coef(pd_model)
pd_logistic <- params[1] / (1 + exp(-params[2] * (0:36 - params[3])))
puma$pd_reconstruction[2:26] <- pd_y
puma$pd_logistic <- pd_logistic


puma$pd_odfw_est[8:35] <- c(
  926,  1058, 1187, 1301, 1393, 1436, 1502, 1590, 1570, 1570, 1592, 1646, 1640, 
  1599, 1592, 1596, 1578, 1541, 1532, 1640, 1703, 1724, 1748, 1760, 1800, 1807, 
  1849, 1910
)
puma$pd_mortalities[1:33] <- c(
  76,   64,  72,  85,  81,  98,  86,  96,  36,  46,  64,  91, 140, 136, 125, 
  142, 149, 168, 140, 165, 177, 171, 158, 162, 169, 164, 135,  94, 110, 113, 
  141, 113, 125
)

pd_mortalities <- tibble(
  year = 1987:2019,
  n = c(76,   64,  72,  85,  81,  98,  86,  96,  36,  46,  64,  91, 140, 136, 
        125, 142, 149, 168, 140, 165, 177, 171, 158, 162, 169, 164, 135,  94, 
        110, 113, 141, 
        113, 125),
  method = "mortalities"
)

all_sources <- puma %>%
  filter(
    !is.na(pd_reconstruction) & 
    !is.na(pd_odfw_est) & 
    !is.na(pd_mortalities)
  ) %>%
  pull(yr)
asi <- which(puma$yr %in% all_sources)
puma <- puma %>%
  mutate(
    pd_recon = scl_rb(pd_reconstruction, c(asi[1], asi[length(asi)])),
    pd_morts = scl_rb(pd_mortalities, c(asi[1], asi[length(asi)])),
    pd_odfwe = scl_rb(pd_odfw_est, c(asi[1], asi[length(asi)])),
    pd_logis = scl_rb(pd_logistic, c(asi[1], asi[length(asi)]))
  ) %>%
  select(yr, pd_recon, pd_morts, pd_odfwe, pd_logis)

# puma %>%
#   pivot_longer(cols = 2:5, names_to = "Source", values_to = "Density") %>%
#   filter(!is.na(Density)) %>%
#   ggplot(aes(x = yr, y = Density, color = Source)) +
#   geom_line()

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
  summarise(ndvi = mean(ndvi)) %>%
  mutate(ndvi = scl(ndvi))

#Prism==========================================================================

prism <- read_csv("data//climate//prism_data_starkey_centroid.csv") %>%
  set_names("dt", "mm_precip", "t_min", "t_mean", "t_max") %>%
  mutate(dt = mdy(dt)) %>%
  mutate(yr = year(dt), mn = month(dt)) %>%
  filter(mn %in% c(5,6,7,8,9)) %>%
  group_by(yr) %>%
  summarise(
    mm_precip = sum(mm_precip),
    t_mean = mean(t_mean)
  ) %>%
  ungroup() %>%
  mutate(precip = scl(mm_precip)) %>%
  mutate(temp = scl(t_mean)) %>%
  select(yr, precip, temp)

# with(prism, plot(t_max ~ t_mean))

#PDI============================================================================

pdi <- read_csv("data//climate//pdi.csv") %>%
  filter(yr >= 1987 & yr < 2024) %>%
  mutate(growing = (may + june + july + august + september) / 5) %>%
  select(yr, growing, september) %>%
  mutate(pdi_september = scl(september)) %>%
  mutate(pdi_growing = scl(growing))

#Clean-up=======================================================================

covariates <- puma %>%
  full_join(pdi) %>%
  full_join(prism) %>%
  full_join(ndvi)

saveRDS(covariates, "data//scaled_covariates.rds")
rm(list = ls())