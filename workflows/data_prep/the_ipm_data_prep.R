# IPM Data Prep
# Kenneth Loonam
# May 2024
# This is the "final" and full version of the IPM data prep for re-submission

#Environment====================================================================

n_yrs <- 34
yr_a <- 1988
yr_z <- 2021
ipm_data <- list()

#Packages=======================================================================

require(dplyr); require(lubridate); require(purrr); require(tidyr)

#CJS data prep==================================================================

load("data//elk_data.Rdata")

y <- elk_data$cap_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()

l <- elk_data$hnt_tib %>%
  arrange(id) %>%
  mutate('2021' = 1) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() %>%
  apply(., 1, function(x) min(which(x != 0)))

w <- elk_data$hnt_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix() 

y <- w + y # harvests count as observed alive that year (did not die naturally)

male <- elk_data$sex_tib %>%
  arrange(id) %>%
  mutate(male = as.numeric(Sex == "M")) %>%
  pull(male)

calf <- elk_data$age_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  as.matrix()
# had a replace_na(0) line appended

herd <- elk_data$hrd_tib %>%
  arrange(id) %>%
  select(as.character(start_year:end_year)) %>%
  mutate_all(funs(case_when(
    . == "main" ~ 0,
    T ~ 1
  ))) %>%
  as.matrix()

gone_elk <- apply(herd, 1, sum) == ncol(herd)
unseen_elk <- apply(y, 1, sum) == 0
bad_elk <- which((gone_elk + unseen_elk) != 0)
f <- apply(y, 1, function(x) min(which(x != 0)))
weird_elk <- f == l


