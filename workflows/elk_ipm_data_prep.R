# IPM data prep -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

# User Specified
n_yrs <- 33
cjs_yrs <- c(min_yr = 1, max_yr = 31)
ratio_yrs <- c(min_yr = 2, max_yr = 31)
first_year <- 1988

cjs_file <- "rand_effect.Rdata"
ratio_file <- "results//recruitment_years_2to31.Rdata"
abundance_file <- "data//elk_abundance_ipm_ready.csv"

# Derived


#Environment====================================================================

require(tidyverse); require(R2jags)

#CJS_Survival===================================================================

load(cjs_file)
cjs_raw <- rslt$BUGSoutput$summary
rm(rslt)

keep_rows <- c(
  cjs_yrs[1]:cjs_yrs[2], 
  cjs_yrs[1]:cjs_yrs[2] + n_yrs, 
  cjs_yrs[1]:cjs_yrs[2] + 2*n_yrs
)

cjs_data <- cjs_raw[grep("phi", dimnames(cjs_raw)[[1]])[keep_rows], ]

s_cjs <- matrix(NA, nrow = nrow(cjs_data), ncol = 5)
s_cjs[,1] <- rep(cjs_yrs[1]:cjs_yrs[2], 3)
s_cjs[,2] <- c(rep(3, nrow(cjs_data)/3*2), 
               rep(2, nrow(cjs_data)/3))
s_cjs[,3] <- c(rep(1, nrow(cjs_data)/3), 
               rep(2, nrow(cjs_data)/3), 
               rep(1, nrow(cjs_data)/3))
s_cjs[,4] <- cjs_data[,1]
s_cjs[,5] <- 1/cjs_data[,2]^2

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
n_sight_ca[,1] <- n_sight_ca[,1] - first_year

n_sight_af <- n_raw %>%
  filter(age == 2 & sex == 1) %>%
  as.matrix(.)
n_sight_af[,1] <- n_sight_af[,1] - first_year

n_sight_am <- n_raw %>%
  filter(age == 2 & sex == 2) %>%
  as.matrix(.)
n_sight_am[,1] <- n_sight_am[,1] - first_year

#Minimum_known_alive============================================================

load("data//elk_minimum_count_data.Rdata")
mindat <- as.matrix(mindat)

#Harvested======================================================================

load("data//elk_harvest_data.Rdata")

#Pasture_changes================================================================

load("data//elk_net_pasture_movement_data.Rdata")
r_df <- movdat %>%
  mutate(mean = mean - hntdat$mean)

rem_dat <- array(NA, dim = c(3,2,34))
for(i in 1:nrow(r_df)){
  rem_dat[r_df$age[i], r_df$sex[i], r_df$year[i]] <- r_df$mean[i]
}
#Combine========================================================================

ipm_data <- list(
  s_cjs = s_cjs,
  r_ratio = r_ratio,
  n_sight_ca = n_sight_ca,
  n_sight_af = n_sight_af,
  n_sight_am = n_sight_am,
  min_n_live = mindat,
  net_remove = rem_dat
)

save(ipm_data, file = "data//elk_ipm_data.Rdata")
