# Pregnancy data prep from full animal handling
# Kenneth Loonam
# June 2023

#Environment====================================================================
# Variables
file_name <- "data//summary_handling_records.csv"
save_file <- "data//pregnancy_data.csv"

# Packages
require(tidyverse); require(lubridate)

#Workflow=======================================================================

rd <- read_csv(file_name, guess_max = 100000) %>%
  filter(Species == "E") %>%
  filter(Sex == "F") %>%
  filter(!is.na(EventDate)) %>%
  mutate(e_date = mdy(EventDate)) %>%
  filter(month(e_date) %in% c(10:12)) %>%
  # filter(WinterFeedgroundDir == "Incoming") %>%
  # filter(!is.na(CalculatedAge)) %>%
  mutate(age_class = case_when(
    CalculatedAge < 3 ~ "young",
    CalculatedAge > 13 ~ "old",
    is.na(CalculatedAge) ~ "prime",
    T ~ "prime"
  )) %>%
  # filter(Lactating %in% c("T", "F")) %>%
  filter(!is.na(Pregnant)) %>%
  filter(Herd == "MAINS")

preg_rates <- rd %>% group_by(HandlingYr, Lactating, age_class) %>%
  summarise(
    pregnancy_rate = mean(Pregnant, na.rm = T),
    pregnancy_sd = sd(Pregnant, na.rm = T),
    n_observations = n()
  ) %>%
  filter(Lactating == "F" | Lactating == "T") %>%
  mutate(lactating = Lactating == "T") %>%
  ungroup() %>%
  separate(HandlingYr, c("yr", "yr2")) %>%
  mutate(id = 1:nrow(.)) %>%
  select(
    id, yr, lactating, pregnancy_rate, pregnancy_sd, n_observations, age_class
    ) 

write.csv(preg_rates, "data//pregnancy_rates.csv")
