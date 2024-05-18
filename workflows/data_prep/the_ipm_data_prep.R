# IPM Data Prep
# Kenneth Loonam
# May 2024
# This is the "final" and full version of the IPM data prep for re-submission

#Environment====================================================================

n_yrs <- 34
yr_a <- 1988
yr_z <- 2021
dtf <- list()

#Packages=======================================================================

require(dplyr); require(purrr); require(tidyr); require(readr)

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


dtf$y <- y; rm(y)
dtf$z <- z; rm(z)
dtf$male <- male; rm(male); rm(male_id)
dtf$female <- female; rm(female)
dtf$calf <- calf; rm(calf)
dtf$herd <- herd; rm(herd)
dtf$l <- l; rm(l)
dtf$f <- f; rm(f)
dtf$nind <- nrow(dtf$y)
rm(gone_elk); rm(unseen_elk); rm(weird_elk); rm(rm_elk); rm(l_k); rm(i); rm(w)
rm(elk_data)

#Ratio data=====================================================================

dtf$r_dt <- read_csv("data//min_n_handle_summaries.csv") %>%
  mutate(yr = year - yr_a + 1) %>%
  select(yr, ca, fe_ad) %>%
  as.matrix() 

#Rebuild from previous IPM data=================================================

load("data//elk_ipm_data_21apr2023.Rdata")

dtf$NC_est <- ipm_data$n_sight_ca %>%
  as_tibble() %>%
  mutate(sd = sqrt(1/tau)) %>%
  select(year, mean, sd) %>%
  as.matrix()
dtf$NF_est <- ipm_data$n_sight_af %>%
  as_tibble() %>%
  mutate(sd = sqrt(1/tau)) %>%
  select(year, mean, sd) %>%
  as.matrix()
dtf$NM_est <- ipm_data$n_sight_am %>%
  as_tibble() %>%
  mutate(sd = sqrt(1/tau)) %>%
  select(year, mean, sd) %>%
  as.matrix()
dtf$NCman <- colSums(ipm_data$n_ca_add - abs(ipm_data$n_ca_rem))
dtf$NFman <- ipm_data$n_ad_add[1,] - abs(ipm_data$n_ad_rem[1,])
dtf$NMman <- ipm_data$n_ad_add[2,] - abs(ipm_data$n_ad_rem[2,])
dtf$n_year <- n_yrs
dtf$nNC <- nrow(dtf$NC_est)
dtf$nNF <- nrow(dtf$NF_est)
dtf$nNM <- nrow(dtf$NM_est)
dtf$nR <- nrow(dtf$r_dt)
dtf$NFhar <- apply(ipm_data$n_hnt[,1,], 2, sum)
dtf$NMhar <- apply(ipm_data$n_hnt[,2,], 2, sum)

alt_min_n <- read_csv("data//min_n_handle_summaries.csv")
NFmin <- ipm_data$min_ad[1,]
for(i in 2:length(NFmin)){
  if(NFmin[i] < alt_min_n$fe_ad[i-1]){
    NFmin[i] <- alt_min_n$fe_ad[i-1]
  }
}
dtf$NFmin <- NFmin; rm(NFmin)
NMmin <- ipm_data$min_ad[2,]
for(i in 2:length(NMmin)){
  if(NMmin[i] < alt_min_n$ma_ad[i-1]){
    NMmin[i] <- alt_min_n$ma_ad[i-1]
  }
}
dtf$NMmin <- NMmin; rm(NMmin)
NCmin <- ipm_data$min_ca
for(i in 2:length(NCmin)){
  if(NCmin[i] < alt_min_n$ca[i-1]){
    NCmin[i] <- alt_min_n$ca[i-1]
  }
}
dtf$NCmin <- NCmin; rm(NCmin); rm(alt_min_n); rm(i)
dtf$NF_ct <- colSums(dtf$y * dtf$female * (1 - dtf$herd), na.rm = T)
dtf$NM_ct <- colSums(dtf$y * dtf$male   * (1 - dtf$herd), na.rm = T)

dtf$nc1e <- ipm_data$est_n1[1,] %>% sum()
dtf$nf1e <- ipm_data$est_n1[2:3,1] %>% sum()
dtf$nm1e <- ipm_data$est_n1[2:3,2] %>% sum()
dtf$nc1e_min <- ipm_data$min_n1[1,] %>% sum()
dtf$nf1e_min <- ipm_data$min_n1[2:3,1] %>% sum()
dtf$nm1e_min <- ipm_data$min_n1[2:3,2] %>% sum()

#Weather covariates=============================================================

scale_clean <- function(x){
  out <- (x - mean(x, na.rm = T))/sd(x, na.rm = T)
  return(out)
}

dtf$sep_pdi <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  filter(year >= yr_a & year <= yr_z) %>%
  pull(september) %>%
  scale_clean()
dtf$may_to_sep_pdi <- read_csv("data//climate//pdi_3508_ne_or.csv") %>%
  filter(year >= yr_a & year <= yr_z) %>%
  select(year, may, june, july, august, september) %>%
  pivot_longer(cols = 2:6) %>%
  group_by(year) %>%
  summarize(pdi = mean(value)) %>%
  pull(pdi) %>%
  scale_clean()

dtf$ndvi_avhrr <- read_csv("data//climate//ndvi_avhrr.csv") %>% 
  filter(yr >= yr_a & yr <= yr_z) %>%
  filter(mn >=5 & mn <= 9) %>%
  group_by(yr) %>%
  summarize(ndvi = mean(NDVI)) %>%
  mutate(ndvi = case_when(
    yr < 2004 ~ ndvi,
    T ~ NA
  )) %>%
  pull(ndvi) %>%
  scale_clean()
dtf$ndvi_modis <- read_csv("data//climate//ndvi_modis.csv") %>% 
  filter(yr >= yr_a & yr <= yr_z) %>%
  filter(mn >=5 & mn <= 9) %>%
  group_by(yr) %>%
  summarize(ndvi = mean(NDVI)) %>%
  pull(ndvi) %>%
  scale_clean()
dtf$ndvi_modis <- c(rep(NA, 12), dtf$ndvi_modis)

dtf$summer_temp <- ipm_data$summer_temp %>% scale_clean()
dtf$summer_precip <- ipm_data$summer_precip %>% scale_clean()

#Community covariates===========================================================

dtf$elk_density <- ipm_data$af_density
dtf$puma_composit <- ipm_data$cougar_density

rm(ipm_data)


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
  add_row(year = 1988:2021, method = "logistic_regression", n = dtf$puma_composit)
# require(ggplot2)
# ggplot(full_puma_data, aes(x = year, y = n, color = method)) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Puma estimates") +
#   scale_color_discrete(
#     name = "Source",
#     labels = c("Souce mean", "Logistic fit", "WMU Mortalities", 
#                "Regional Estimate", "WMU Reconstruction")
#   ) +
#   theme_classic() +
#   xlab("Year") +
#   ylab("Scaled estimate")
  
dtf$puma_mean <- full_puma_data %>% 
  filter(method == "group_mean") %>%
  pull(n) %>%
  scale_clean()

saveRDS(dtf, file = "data//the_ipm_data.rds")
rm(list = ls())
