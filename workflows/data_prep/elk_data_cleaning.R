# Starkey Capture and handling data cleaning
# Kenneth Loonam
# July 2020
# Prep data for components of the ipm and save as r objects

#Variables======================================================================

target_species <- species
start_year <- years[1]
end_year <- years[length(years)]
cap_window <- c(11, 3) # captures from november to march count for cjs

#Environment====================================================================

require(tidyverse); require(clock); require(lubridate)

#Initial_data_prep==============================================================

raw_data <- read_csv("data//handling.csv", guess_max = 26877) %>%
  filter(Species == target_species) %>%
  filter(!is.na(Sex)) %>% # excludes 2 elk with no gender listed
  mutate(event_date = mdy(EventDate)) %>%
  mutate(death_date = mdy(DeathDate)) %>%
  mutate(birth_date = mdy(BirthDate)) %>%
  mutate(entry_date = mdy(AnimalHandling_CreateDate)) %>%
  mutate(id = AnimalHandling_AnimalId) %>%
  mutate(event_date = case_when(
    event_date > entry_date ~ entry_date,
    T ~ event_date # events aren't recorded before they happen
  )) 

# rewrite these as a list and the mutate functions as lapply calls?
mains_pastures <- c("BALLY", "DOUGP", "HFMOO", "MAINS")
neeas_pastures <- c("NEEAS")
newes_pastures <- c("NEWES")
nestd_pastures <- c("NESTD")
outsd_pastures <- c("ELKHO", "LADDM", "NCASI", "OUTSD")
campp_pastures <- c("CAMPP")
handl_pastures <- c("BARNP", "BCALY", "BEARP", "BULLP", "CUHNA", "FEEDG", 
  "HANDL", "MCALY", "MDCRE", "MDCRW", "SOUTH", "SPARE", "UPPER", "WINGP")

full_data <- raw_data %>%
  mutate(capture_herd = case_when(
    CapturePasture %in% mains_pastures ~ "main",
    CapturePasture %in% neeas_pastures ~ "ne_e",
    CapturePasture %in% nestd_pastures ~ "ne_o",
    CapturePasture %in% newes_pastures ~ "ne_w",
    CapturePasture %in% campp_pastures ~ "camp",
    CapturePasture %in% outsd_pastures ~ "otsd",
    CapturePasture %in% handl_pastures ~ "hand",
  )) %>%
  mutate(release_herd = case_when(
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
    !(is.na(herd) & is.na(capture_herd) & is.na(release_herd))
  ) %>%
  mutate(population = case_when(
    (!is.na(release_herd) & release_herd != "hand") ~ release_herd,
    (!is.na(capture_herd) & capture_herd != "hand") ~ capture_herd,
    (!is.na(herd) & herd != "hand") ~ herd,
    !is.na(capture_herd) ~ capture_herd,
    !is.na(release_herd) ~ release_herd,
    !is.na(herd) ~ herd
  )) %>%
  mutate(next_bio_year = case_when(
    month(event_date) > 10 ~ year(event_date) + 1,
    T ~ year(event_date)
  ))

x <- full_data %>%
  filter(population == "main") %>%
  filter(CaptureMethod == "HUNTED")

# Population assignment rules:
  # Priority: release_herd, capture_herd, herd
  # NA -> next step in priority
  # "hand" -> next step in priority
    # "hand" only the population if no other herd is given
  # Event interpreted as: individual switched to population "new"
    # If population "new" == "old", no switch occurred, population maintained

#Herd_matrix====================================================================

build_herd_history <- function(targ_id, full_data, start_year, end_year){
  # Assigns an individual to a heard on each Nov 1 after its first observation
  
  herd_id <- tibble(
    year = start_year:(end_year + 1),
    herd = NA
  )
  
  ind_data <- full_data %>%
    filter(id == targ_id) %>%
    arrange(event_date)
  
  for(i in 1:nrow(ind_data)){
    herd_id$herd[herd_id$year == ind_data$next_bio_year[i]] <- ind_data$population[i]
  }
  
  first_obs <- which.min(is.na(herd_id$herd))
  if(first_obs > 1){
    herd_id$herd[first_obs - 1] <- herd_id$herd[first_obs]
  }
  
  for(i in first_obs:nrow(herd_id)){
    if(is.na(herd_id$herd[i])){
      herd_id$herd[i] <- herd_id$herd[i-1]
    }
  }
  out <- c(targ_id, herd_id$herd)
  names(out) <- c("id", herd_id$year)
  return(out)
}

herd_tib <- bind_rows(map(
  .x = unique(full_data$id), 
  .f = build_herd_history, 
  full_data = full_data,
  start_year = start_year, 
  end_year = end_year
))

#Known_alive====================================================================

# alive, dead, unknown

f_l_data <- full_data %>%
  filter(
    CaptureMethod %in% c(
      "TRAPPED", "HUNTED", "HANDLD", "MRT-TRAP", "NETTED", "DARTED", "CLOVER", 
      "MOVEDA",  "OBSRVD", "SEARCH" 
      )
    ) %>%
  group_by(id) %>%
  summarize(
    first_capture = min(event_date),
    last_capture = max(event_date)
  ) %>%
  mutate(first_alive = case_when(
    month(first_capture) <= cap_window[2] ~ year(first_capture) - 1,
    T ~ year(first_capture)
  )) %>%
  mutate(last_alive = case_when(
    month(last_capture) >= cap_window[1] ~ year(last_capture),
    T ~ year(last_capture) - 1
  )) %>%
  ungroup()

alive_tib <- matrix(NA, nrow = nrow(herd_tib), ncol = ncol(herd_tib))

for(i in 1:nrow(herd_tib)){
  
  alive_tib[i,1] <- herd_tib$id[i]
  
  if(herd_tib$id[i] %in% f_l_data$id){
    first_year <- f_l_data[f_l_data$id == herd_tib$id[i],]$first_alive
    if(first_year < start_year){first_year <- start_year} 
    last_year  <- f_l_data[f_l_data$id == herd_tib$id[i],]$last_alive
    if(last_year < start_year){last_year <- start_year} 
    alive_tib[i,(first_year:last_year) - start_year + 2] <- 1
  }
}

alive_tib <- alive_tib %>% as_tibble()
names(alive_tib) <- c("id", as.character(start_year:(end_year+1)))
alive_tib <- alive_tib %>%
  arrange(id)

death_dates <- full_data %>%
  group_by(id) %>%
  summarise(death_date = max(death_date)) %>%
  ungroup() %>%
  arrange(id) %>%
  mutate(dead_year = case_when(
    month(death_date) >= cap_window[1] ~ year(death_date) + 1,
    T ~ year(death_date)
  ))

for(i in 1:nrow(alive_tib)){
  if(!is.na(death_dates$dead_year[i])){
    assertthat::assert_that(death_dates$id[i] == alive_tib$id[i])
    death_column <- death_dates$dead_year[i] - start_year + 2
    if(death_column < 2){death_column <- 2}
    alive_tib[i,death_column:ncol(alive_tib)] <- "0"
  }
}

alive_tib <- alive_tib %>%
  mutate_at(vars(-id), funs(as.numeric))


#Sex============================================================================

sex_tib <- full_data %>%
  distinct(id, Sex) %>%
  arrange(id)

#Age============================================================================

get_birth_date <- function(x, id_x){
  # calculates birthday for individuals
  # only concerned with getting down to adult/calf distinctions
  
  tmp <- x %>%
    filter(id == id_x) %>%
    arrange(event_date) %>%
    select(birth_date, Class, AgeYr, AgeMo, AgeDay, event_date, id)
  
  if(any(!is.na(tmp$birth_date))){ # if birthday is given, use that
    out_date <- min(tmp$birth_date)
  }else{
    if(all(is.na(tmp$Class))){ # if there is no age info
      birth_year <- year(min(tmp$event_date)) - 3 # just assign an arbitrary adult age
      out_date <- ymd(paste0(birth_year, "-06-01"))
    }else{
      if(any(tmp$Class == "CA" | tmp$Class == "NC", na.rm = T)){ # if its captured as a calf
        birth_year <- year(min(tmp$event_date)) - (month(min(tmp$event_date)) < 5)
        out_date <- ymd(paste0(birth_year, "-06-01")) # use most recent June 1st
      }else{
        if(any(tmp$Class == "AD" | tmp$Class == "AG", na.rm = T)){ # if captured adult
          birth_year <- year(min(tmp$event_date)) - (month(min(tmp$event_date)) < 5) - 1
          out_date <- ymd(paste0(birth_year, "-06-01")) # 2nd most recent June 1st
        }else{ # If they're older than that it won't matter. They haven't been observed before,
          out_date <- NA # so they won't be used in any of the analyses for those years.
        }
      }
    }
  }
  
  out <- tibble(id = id_x, birth_date = out_date)
}


age_data <- full_data %>% 
  map(.x = unique(.$id), .f = get_birth_date, x = .) %>% # rewrite this function
  bind_rows() %>%
  arrange(id) %>%
  mutate(birth_date = case_when(
    id == "151203E01" ~ mdy("6-1-2015"),
    id != "151203E01" ~ .$birth_date)) %>%
  mutate(calf_session = year(birth_date) - start_year + 1)
age_data <- age_data %>%
  mutate(calf_session = case_when(
    calf_session < 0  ~ 0,
    calf_session >= 0 ~ calf_session
  ))

age_tib <- matrix(NA, nrow = nrow(age_data), ncol = ncol(herd_tib) - 1)
for(i in 1:nrow(age_data)){
  age_tib[i, age_data$calf_session[i]] <- 1
}
for(i in 1:nrow(age_tib)){
  if(all(is.na(age_tib[i,]))){
    age_tib[i,] <- 0
  }else{
    age_tib[i, (1 + which.min(age_tib[i,] == 1)):ncol(age_tib)] <- 0
  }
}

age_tib <- age_data %>%
  select(id) %>%
  bind_cols(as_tibble(age_tib))
names(age_tib) <- c("id", as.character(start_year:(end_year + 1)))

#Capture_history================================================================

cap_data <- full_data %>%
  filter(CaptureMethod %in% c(
    "TRAPPED", "HANDLD", "NETTED", "DARTED", "CLOVER", "MOVEDA",  "OBSRVD")) %>%
  filter(event_date >= ymd(paste(start_year, cap_window[1], "1"))) %>%
  mutate(cap_month = month(event_date)) %>%
  filter(cap_month >= cap_window[1] | cap_month <= cap_window[2]) %>%
  mutate(session = case_when(
    cap_month >= cap_window[1] ~ year(event_date) - start_year + 1,
    cap_month <= cap_window[2] ~ year(event_date) - start_year
  )) %>%
  arrange(id, session)

cap_tib <- matrix(
  0, 
  nrow = length(unique(cap_data$id)), 
  ncol = ncol(alive_tib) - 1) %>%
  as_tibble() %>%
  bind_cols(unique(cap_data$id), .) 
names(cap_tib) <- c("id", as.character(start_year:(end_year + 1)))

for(i in 1:nrow(cap_data)){
  cap_tib[cap_tib$id == cap_data$id[i], cap_data$session[i] + 1] <- 1
}

#Hunt_history===================================================================

hunt_data <- full_data %>%
  filter(CaptureMethod == "HUNTED") %>%
  mutate(cap_month = month(event_date)) %>%
  mutate(session = case_when(
    cap_month >= cap_window[1] ~ year(event_date) - start_year + 2,
    T ~ year(event_date) - start_year + 1
  )) %>%
  arrange(id, session)

hunt_tib <- matrix(
  0, 
  nrow = length(unique(alive_tib$id)), 
  ncol = ncol(alive_tib) - 1) %>%
  as_tibble() %>%
  bind_cols(unique(alive_tib$id), .) 
names(hunt_tib) <- c("id", as.character(start_year:(end_year + 1)))

for(i in 1:nrow(hunt_data)){
  hunt_tib[hunt_tib$id == hunt_data$id[i], hunt_data$session[i] + 1] <- 1
}

#Clean_up=======================================================================

elk_data <- list(
  cap_tib = cap_tib,
  age_tib = age_tib,
  sex_tib = sex_tib,
  liv_tib = alive_tib,
  hrd_tib = herd_tib,
  hnt_tib = hunt_tib,
  explainer = "cap is capture history for cjs,
  age is year that id was a calf,
  sex is sex of id,
  liv is 1 for individual known alive, 0 for dead, NA for unknown,
  hrd is which herd id belonged to on 11-1 of each year,
  hnt is a 1 for the year an elk was removed from the population (hunting)"
)

save(elk_data, file = elk_data_location)

rm(list = ls())
