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
  filter(value != 0) %>%
  mutate(year = year - 1987)
mindat <- matrix(1, nrow = nrow(min_dat), ncol = 5)
mindat[,1] <- min_dat$year
mindat[which(min_dat$name == "ca"),2] <- 1
mindat[which(min_dat$name != "ca"),2] <- 3


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

save(ipm_data, file = "data//elk_ipm_data_survival_includes_harvest.Rdata")
