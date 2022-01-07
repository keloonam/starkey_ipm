# IPM data prep -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# User Specified
n_yrs <- 33
cjs_yrs <- c(min_yr = 1, max_yr = 31)
ratio_yrs <- c(min_yr = 2, max_yr = 31)
first_year <- 1988

cjs_file <- "results//survival//cjs_rslt_03jan2022.Rdata"
ratio_file <- "results//recruitment_years_2to31.Rdata"
abundance_file <- "data//elk_abundance_ipm_ready.csv"

# Survival estimates to remove from CJS
# Parameters either did not converge or 
# are from years with unrecorded id's removed
cjs_removals <- c(
  "af[1]",
  "af[3]",
  "af[32]",
  "af[33]",
  "am[1]",
  "am[2]",
  "am[3]",
  "am[4]",
  "am[5]",
  "am[6]",
  "am[6]",
  "am[8]",
  "am[10]",
  "am[12]",
  "am[18]",
  "am[31]",
  "am[32]",
  "am[33]",
  "ca[1]",
  "ca[4]",
  "ca[17]",
  "ca[32]",
  "ca[33]",
  "af[29]",
  "af[11]",
  "af[13]",
  "af[14]"
)


#Environment====================================================================

require(tidyverse)

#CJS_Survival===================================================================

load(cjs_file)
cjs_raw <- cjs_rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(starts_with("survival_")) %>%
  select(-contains(cjs_removals))
rm(cjs_rslt)

# first year is survival to 1989
cjs_years <- as.numeric(str_extract_all(names(cjs_raw), "[0-9]+"))
cjs_males <- grep("am", names(cjs_raw))
cjs_calfs <- grep("ca", names(cjs_raw))
cjs_means <- unlist(map(cjs_raw, mean))
cjs_sd    <- unlist(map(cjs_raw, sd))

cjs_dat <- matrix(1, nrow = length(cjs_years), ncol = 5)
cjs_dat[,1] <- cjs_years + 1
cjs_dat[,2] <- 3
cjs_dat[cjs_calfs,2] <- 2
cjs_dat[cjs_males,3] <- 2
cjs_dat[,4] <- cjs_means
cjs_dat[,5] <- 1/(cjs_sd^2)
cjs_dat <- cjs_dat[-69,]
#Ratio_Recruitment==============================================================

load(ratio_file)
ratio_raw <- rslt$BUGSoutput$summary
rm(rslt)

ratio_data <- ratio_raw[grep("R", dimnames(ratio_raw)[[1]]),]

r_ratio <- matrix(NA, nrow = nrow(ratio_data), ncol = 5)
r_ratio[,1] <- ratio_yrs[1]:ratio_yrs[2]
r_ratio[,4] <- ratio_data[,1]
r_ratio[,5] <- 1/ratio_data[,2]^2

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

#Minimum_known_alive============================================================

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

#Harvested======================================================================

load("data//elk_harvest_data.Rdata")
n_hnt <- array(data = 0, dim = c(4,2,n_yrs))
hntdat <- as.matrix(hntdat)
for(i in 1:(nrow(hntdat) - 6)){
  n_hnt[hntdat[i,2] + 1, hntdat[i,3], hntdat[i,1]] <- hntdat[i,4]
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

n_ad_rem <- n_rem[2,,] + n_rem[3,,]
n_ca_rem <- n_rem[1,,]
n_ad_add <- n_add[2,,] + n_rem[3,,]
n_ca_add <- n_add[1,,]

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
  n_hnt = n_hnt
)

save(ipm_data, file = "data//elk_ipm_data_07jan2022.Rdata")
