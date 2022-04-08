# Number Harvested -- Starkey Elk
# Kenneth Loonam
# April 2021

#Variables======================================================================

target_herd <- target_herd
start_year <- years[1]
end_year <- years[length(years)]
n_yr <- length(years)

#Environment====================================================================

require(tidyverse)
load("data//elk_data.Rdata")

#Data_prep======================================================================

# Get ages in jags format (1 = calf, 2 = yearling, 3 = subadult, 4 = adult)
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
hnt_tib <- elk_data$hnt_tib #%>%
  #mutate(across(where(is.numeric), ~ replace_na(., 0)))
hrd_tib <- elk_data$hrd_tib %>%
  mutate(across(where(is.character), ~ replace_na(., target_herd)))


# Prepare empty minimum count data tibble
hntdat <- tibble(
  year = sort(rep((start_year:end_year) - start_year + 1, 8)),
  age = rep(c(1,2,3,4), n_yr*2),
  sex = rep(c(1,1,1,1,2,2,2,2), n_yr),
  mean = 0
)

# Fill the tibble
for(j in 1:nrow(hntdat)){
  yr <-  as.numeric(hntdat[j,1])
  age <- as.numeric(hntdat[j,2])
  sex <- as.numeric(hntdat[j,3])
  
  hntdat[j,4] <- sum((age_tib[,yr+1] == age) * (sex_tib$sex == sex) * (hrd_tib[,yr+1] == target_herd) * hnt_tib[,yr+1])
}

# bob <- hntdat %>%
#   group_by(year) %>%
#   summarize(n = sum(mean)) %>%
#   ungroup()
# plot(bob)
# hist(bob$n)

save(hntdat, file = "data//elk_harvest_data.Rdata")
