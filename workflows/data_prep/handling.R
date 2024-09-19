# Clean the capture and handling data to prep ipm data
# Kenneth Loonam
# August 2024

#Environment====================================================================
# Define Pastures
mains_pastures <- c("BALLY", "DOUGP", "HFMOO", "MAINS")
neeas_pastures <- c("NEEAS")
newes_pastures <- c("NEWES")
nestd_pastures <- c("NESTD")
outsd_pastures <- c("ELKHO", "LADDM", "NCASI", "OUTSD")
campp_pastures <- c("CAMPP")
handl_pastures <- c("BARNP", "BCALY", "BEARP", "BULLP", "CUHNA", "FEEDG", 
                    "HANDL", "MCALY", "MDCRE", "MDCRW", "SOUTH", "SPARE", 
                    "UPPER", "WINGP")

#Prep_main_data_tibble==========================================================
rd <- read_csv(ah_fp, guess_max = 100000) %>%
  filter(Species == "E") %>% # filter to elk
  filter(Sex == "F" | Sex == "M") %>% # exclude entries without known sex
  mutate(event_dt = mdy(EventDate)) %>% 
  mutate(death_dt = mdy(DeathDate)) %>%
  mutate(birth_dt = mdy(BirthDate)) %>%
  mutate(entry_dt = mdy(CreateDate)) %>%
  mutate(id = AnimalHandling_AnimalId) %>%
  mutate(event_dt = case_when(
    event_dt > entry_dt ~ entry_dt,
    T ~ event_dt # events aren't recorded before they happen
  )) %>%
  mutate(cap_herd = case_when(
    CapturePasture %in% mains_pastures ~ "main",
    CapturePasture %in% neeas_pastures ~ "ne_e",
    CapturePasture %in% nestd_pastures ~ "ne_o",
    CapturePasture %in% newes_pastures ~ "ne_w",
    CapturePasture %in% campp_pastures ~ "camp",
    CapturePasture %in% outsd_pastures ~ "otsd",
    CapturePasture %in% handl_pastures ~ "hand",
  )) %>%
  mutate(rel_herd = case_when(
    ReleasePasture %in% mains_pastures ~ "main",
    ReleasePasture %in% neeas_pastures ~ "ne_e",
    ReleasePasture %in% nestd_pastures ~ "ne_o",
    ReleasePasture %in% newes_pastures ~ "ne_w",
    ReleasePasture %in% campp_pastures ~ "camp",
    ReleasePasture %in% outsd_pastures ~ "otsd",
    ReleasePasture %in% handl_pastures ~ "hand",
  )) %>%
  mutate(herd = case_when(
    Herd %in% mains_pastures ~ "main",
    Herd %in% neeas_pastures ~ "ne_e",
    Herd %in% newes_pastures ~ "ne_w",
    Herd %in% nestd_pastures ~ "ne_o",
    Herd %in% outsd_pastures ~ "otsd",
    Herd %in% campp_pastures ~ "camp",
    Herd %in% handl_pastures ~ "hand",
  )) %>%
  filter( # remove individuals with no herd data
    !(is.na(herd) & is.na(cap_herd) & is.na(rel_herd))
  ) %>%
  mutate(pop = case_when( # assign to a population, priority release data
    (!is.na(rel_herd) & rel_herd != "hand") ~ rel_herd,
    (!is.na(cap_herd) & cap_herd != "hand") ~ cap_herd,
    (!is.na(herd)     & herd     != "hand") ~ herd,
    !is.na(rel_herd)                        ~ rel_herd,
    !is.na(cap_herd)                        ~ cap_herd,
    !is.na(herd)                            ~ herd
  )) %>%
  mutate(next_year = case_when( # what is the next Nov1 that will occur?
    month(event_dt) > 10 ~ year(event_dt) + 1,
    T ~ year(event_dt)
  ))

#Prep_herd_info=================================================================

herd_tib <- map( # apply over all ids in full data
  .x = unique(rd$id), 
  .f = build_herd_history, 
  full_data = rd,
  start_year = start_year, 
  end_year = end_year
) %>% bind_rows() %>% arrange(id)

#Prep_known_alive===============================================================

known_alive_tib <- map( # apply over all ids in full data
  .x = unique(rd$id), 
  .f = build_known_alive_tibble, 
  full_data = rd,
  start_yr = start_year, 
  end_yr = end_year
) %>% bind_rows() %>% arrange(id) %>%
  mutate(across(2:ncol(.), as.numeric))

#Prep_sex_tibble================================================================

sex_tib <- rd %>%
  mutate(sex = Sex) %>%
  distinct(id, sex) %>%
  arrange(id)

#Prep_age_tibble================================================================

age_tib <- map( # apply over all ids in full data
  .x = unique(rd$id), 
  .f = build_age_history, 
  full_data = rd,
  start_yr = start_year, 
  end_yr = end_year
) %>% bind_rows() %>% arrange(id) %>%
  mutate(across(2:ncol(.), as.numeric))

#Prep_capture_history===========================================================

capture_tib <- build_capture_history(
  full_data = rd,
  start_yr = start_year,
  end_yr = end_year,
  cap_window = c(11, 3) # captures from november through march count
) %>% arrange(id)

#Prep_harvest_history===========================================================

harvest_tib <- build_harvest_history(
  full_data = rd,
  start_yr = start_year,
  end_yr = end_year
) %>% arrange(id)

#Save_and_clean_up==============================================================

handling_data <- list(
  capture_history = capture_tib,
  herd_assignment = herd_tib,
  annual_age = age_tib,
  known_alive = known_alive_tib,
  sex = sex_tib,
  harvest_year = harvest_tib
)

saveRDS(handling_data, file = "data//capture_handling_data.rds")
rm(list = ls())
