# Starkey elk CJS model data prep workflow
# Kenneth Loonam
# March 2020

# We need:
  # Capture history from feedgrounds only (matrix)
  # True observed state based on all observations (matrix)
  # Vector of sexes
  # Matrix of ages
  # Matrix of herds
  # Vector of first captures
  # nocc and nind

#Variables======================================================================

target_herd <- "MAINS"
model <- "cjs_phiyear_pyear"
start_year <- 1988 # first captures included at end of fall start_year
cap_window <- c(11, 3) # include captures from month ## to month ## inclusive

#Environment====================================================================

require(tidyverse); require(lubridate)
source("r//elk_cjs_functions.R")

#Data===========================================================================

raw_data <- read_csv("data//elk_handling.csv", guess_max = 26877) %>%
  filter(Species == "E") %>%
  filter(!is.na(Sex)) %>% # excludes 2 elk with no gender listed
  mutate(event_date = mdy(EventDate)) %>%
  mutate(death_date = mdy(DeathDate)) %>%
  mutate(birth_date = mdy(BirthDate)) %>%
  mutate(id = AnimalHandling_AnimalId)

# bob <- raw_data %>% # error checking, I think, could comment out??
#   filter(CaptureMethod == "HANDLD") %>%
#   filter(event_date >= ymd(paste(start_year, cap_window[1], "1"))) %>%
#   mutate(cap_month = month(event_date)) %>%
#   filter(cap_month >= cap_window[1] | cap_month <= cap_window[2]) %>%
#   mutate(session = case_when(
#     cap_month >= cap_window[1] ~ year(event_date) - start_year + 1,
#     cap_month <= cap_window[2] ~ year(event_date) - start_year
#   )) %>%
#   arrange(id, session) %>%
#   mutate(id3 = as.numeric(as.factor(id))) %>%
#   filter(id3 == 1441)
  
  

# Record captured individual and session number
cap_data <- raw_data %>%
  filter(CaptureMethod == "HANDLD") %>%
  filter(event_date >= ymd(paste(start_year, cap_window[1], "1"))) %>%
  mutate(cap_month = month(event_date)) %>%
  filter(cap_month >= cap_window[1] | cap_month <= cap_window[2]) %>%
  mutate(session = case_when(
    cap_month >= cap_window[1] ~ year(event_date) - start_year + 1,
    cap_month <= cap_window[2] ~ year(event_date) - start_year
  )) %>%
  arrange(id, session)

# Record sex of captured individuals
sex_data <- raw_data %>%
  distinct(id, Sex) %>%
  filter(id %in% cap_data$id) %>%
  arrange(id) %>%
  mutate(id = as.numeric(as.factor(id)))

# Record birth dates of captured ids
# assumes 185 elk with no data were adults at first capture
# Did not check if those ids made it into final capture list
age_data <- raw_data %>% 
  map(.x = unique(.$id), .f = get_birth_date, x = .) %>% # rewrite this function
  bind_rows() %>%
  filter(id %in% cap_data$id) %>%
  arrange(id) %>%
  mutate(birth_date = case_when(
    id == "151203E01" ~ mdy("6-1-2015"),
    id != "151203E01" ~ .$birth_date)) %>%
  mutate(id = as.numeric(as.factor(id))) %>%
  mutate(calf_session = year(birth_date) - start_year + 1)
age_data <- age_data %>%
  mutate(calf_session = case_when(
    calf_session < 0  ~ 0,
    calf_session >= 0 ~ calf_session
  )) %>%
  filter(calf_session > 0)

# Standardize cap_data id column with other data (1:n_id)
cap_data <- cap_data %>%
  mutate(id = as.numeric(as.factor(id))) %>%
  select(id, session)

#Matrices=======================================================================
# Prep data for JAGS (arrays\vectors\matrices, numeric)

cap_mat <- matrix(0, nrow = nrow(sex_data), ncol = max(cap_data$session))
for(i in 1:nrow(cap_data)){
  cap_mat[cap_data$id[i], cap_data$session[i]] <- 1
}

sex_mat <- as.numeric(sex_data$Sex == "M")

# Everything I have done with ages in the project is hideous >.<
age_mat <- matrix(NA, nrow = nrow(sex_data), ncol = max(cap_data$session))
for(i in 1:nrow(age_data)){
  age_mat[age_data$id[i], age_data$calf_session[i]] <- 1
}
for(i in 1:nrow(age_mat)){
  if(all(is.na(age_mat[i,]))){
    age_mat[i,] <- 0
  }else{
    age_mat[i, (1 + which.min(age_mat[i,] == 1)):ncol(age_mat)] <- 0
  }
}

f <- apply(cap_mat, 1, function(x) min(which(x != 0)))

elk_cjs_data <- list(
  n_occ = ncol(cap_mat),
  n_ind = nrow(cap_mat),
  y = cap_mat,
  f = f,
  male = sex_mat,
  calf = age_mat
)

save(elk_cjs_data, file = "data//elk_cjs_data.Rdata")
