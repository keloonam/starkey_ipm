# Starkey Elk minimum number known alive workflow
# Kenneth Loonam
# June 2020

#Variables======================================================================

targe_herd <- target_herd
start_year <- years[1]
end_year <- years[length(years)]
n_yr <- length(years)

#Environment====================================================================

require(tidyverse)

load("data//elk_data.Rdata")
# raw_data <- read_csv("handling.csv")
fdg_data <- read_csv("data//min_n_handle_summaries.csv")

#Workflow=======================================================================

# Get ages in jags format (1 = calf, 2 = subadult, 3 = adult)
age_tib <- elk_data$age_tib
tmp <- apply(
  X = as.matrix(age_tib[,2:ncol(age_tib)]),
  MARGIN = 1,
  FUN = function(x) min(which(x == 1))
)
tmp[tmp == Inf] <- NA
for(i in 1:nrow(age_tib)){
  if(!is.na(tmp[i])){
    age_tib[i,tmp[i] + 2] <- 2
    age_tib[i,tmp[i] + 3] <- 3
    age_tib[i,(tmp[i] + 4):ncol(age_tib)] <- 4
  }
  if(all(is.na(age_tib[i,2:ncol(age_tib)]))){
    age_tib[i,2:ncol(age_tib)] <- 4
  }
}
age_tib <- age_tib %>%
  mutate(across(where(is.numeric), ~ replace_na(., 4))) %>%
  mutate(across(where(is.numeric), ~ if_else(.==0, 4, .)))

# Get sex in jags format (1 = F, 2 = M)
sex_tib <- elk_data$sex_tib %>%
  mutate(sex = ((Sex == "M") + 1)) %>%
  select(id, sex)

# Pull out other tibbles
liv_tib <- elk_data$liv_tib %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0)))
hrd_tib <- elk_data$hrd_tib %>%
  mutate(across(where(is.character), ~ replace_na(., "no_herd")))

# Prepare empty minimum count data tibble
mindat <- tibble(
  year = sort(rep((start_year:end_year) - start_year + 1, 8)),
  age = rep(c(1,2,3,4), n_yr*2),
  sex = rep(c(1,1,1,1,2,2,2,2), n_yr),
  mean = 0
)

# Fill the tibble
for(i in 1:nrow(mindat)){
  yr <-  as.numeric(mindat[i,1])
  age <- as.numeric(mindat[i,2])
  sex <- as.numeric(mindat[i,3])
  
  mindat[i,4] <- sum((age_tib[,yr+1] == age) * (sex_tib$sex == sex) * liv_tib[,yr+1] * (hrd_tib[,yr+1] == "main"))
}

mindat_af <- mindat %>%
  filter(age > 1) %>%
  filter(sex == 1) %>%
  group_by(year) %>%
  summarise(af = sum(mean)) %>%
  mutate(year = year + 1987)
mindat_am <- mindat %>%
  filter(age > 1) %>%
  filter(sex == 2) %>%
  group_by(year) %>%
  summarise(am = sum(mean)) %>%
  mutate(year = year + 1987)
mindat_ca <- mindat %>%
  filter(age == 1) %>%
  group_by(year) %>%
  summarise(c = sum(mean)) %>%
  mutate(year = year + 1987)

min_counts <- mindat_af %>%
  left_join(mindat_am) %>%
  left_join(mindat_ca) %>%
  left_join(fdg_data) %>%
  replace_na(list(fe_ad = 0, ma_ad = 0, ca = 0)) %>%
  mutate(af = case_when(fe_ad > af ~ fe_ad, T ~ af)) %>%
  mutate(am = case_when(ma_ad > am ~ ma_ad, T ~ am)) %>%
  mutate(c = case_when(ca > c ~ ca, T ~ c)) %>%
  select(1:4) %>%
  rename(ca = c)

save(min_counts, file = "data//elk_minimum_count_data.Rdata")
