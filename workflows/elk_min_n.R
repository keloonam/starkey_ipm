# Starkey Elk minimum number known alive workflow
# Kenneth Loonam
# June 2020

#Variables======================================================================

targe_herd <- "main"
start_year <- 1988
end_year <- 2021

#Environment====================================================================

require(tidyverse)

load("data//elk_data.Rdata")

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
    age_tib[i,(tmp[i] + 3):ncol(age_tib)] <- 3
  }
  if(all(is.na(age_tib[i,2:ncol(age_tib)]))){
    age_tib[i,2:ncol(age_tib)] <- 3
  }
}
age_tib <- age_tib %>%
  mutate(across(where(is.numeric), ~ replace_na(., 3))) %>%
  mutate(across(where(is.numeric), ~ if_else(.==0, 3, .)))

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
  year = sort(rep((start_year:end_year) - start_year + 1, 6)),
  age = rep(c(1,2,3), 34*2),
  sex = rep(c(1,1,1,2,2,2), 34),
  mean = 0
)

# Fill the tibble
for(i in 1:nrow(mindat)){
  yr <-  as.numeric(mindat[i,1])
  age <- as.numeric(mindat[i,2])
  sex <- as.numeric(mindat[i,3])
  
  mindat[i,4] <- sum((age_tib[,yr+1] == age) * (sex_tib$sex == sex) * liv_tib[,yr+1] * (hrd_tib[,yr+1] == "main"))
}

x <- mindat %>%
  group_by(year) %>%
  summarize(n = sum(mean))
plot(x)

save(mindat, file = "data//elk_minimum_count_data.Rdata")
