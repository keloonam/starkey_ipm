# Number of elk moved -- Starkey IPM
# Kenneth Loonam
# August 2021

#Variables======================================================================

target_herd <- "main"
start_year <- 1988
end_year <- 2021

#Environment====================================================================

require(tidyverse)
load("data//elk_data.Rdata")

#Data_management================================================================

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

#Workflow=======================================================================

herd_change <- matrix(0, nrow = nrow(hrd_tib), ncol = ncol(hrd_tib) - 1)

for(i in 1:nrow(hrd_tib)){
  for(j in 3:ncol(hrd_tib)){
    k <- j-1
    if(hrd_tib[i,j] == "main"){
      if(hrd_tib[i,j-1] != "main" & hrd_tib[i,j-1] != "no_herd"){
        herd_change[i,k] <- 1 # added
      }
      if(j < ncol(hrd_tib)){
        if(hrd_tib[i,j+1] != "main" & hrd_tib[i,j+1] != "no_herd"){
          herd_change[i,k+1] <- -1 # removed
        }
      }
    }
  }
}

# Prepare empty net movement data tibble
movdat <- tibble(
  year = sort(rep((start_year:end_year) - start_year + 1, 6)),
  age = rep(c(1,2,3), 34*2),
  sex = rep(c(1,1,1,2,2,2), 34),
  mean = 0
)

# Fill the tibble
for(i in 1:nrow(movdat)){
  yr <-  as.numeric(movdat[i,1])
  age <- as.numeric(movdat[i,2])
  sex <- as.numeric(movdat[i,3])
  
  movdat[i,4] <- sum((age_tib[,yr+1] == age) * (sex_tib$sex == sex) * liv_tib[,yr+1] * herd_change[,yr])
}

x <- movdat %>%
  group_by(year) %>%
  summarize(n = sum(mean))
plot(x)

save(movdat, file = "data//elk_net_pasture_movement_data.Rdata")
