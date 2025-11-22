# Starkey Capture and handling data cleaning
# Kenneth Loonam
# July 2020
# Prep data for components of the ipm and save as r objects

#Variables======================================================================

elk_data_location <- "data//elk_data.Rdata"
target_species <- "E"
start_year <- 1988
end_year <- 2021
cap_window <- c(11, 3) # captures from November to March count for cjs

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
  mutate(next_bio_year = case_when( # what is the next Nov1 that will occur?
    month(event_date) > 10 ~ year(event_date) + 1,
    T ~ year(event_date)
  ))


#Herd_matrix====================================================================
# Population assignment rules:
# Priority: release_herd, capture_herd, herd
# NA -> next step in priority
#  "hand" -> next step in priority
#  "hand" only the population if no other herd is given
# Event interpreted as: individual switched to population "new"
# If population "new" == "old", no switch occurred, population maintained

build_herd_history <- function(targ_id, full_data, start_year, end_year){
  # Assigns an individual to a heard on each Nov 1 after its first observation
  
  herd_id <- tibble( # empty single-individual tibble
    year = start_year:(end_year + 1),
    herd = NA
  )
  
  ind_data <- full_data %>% # get events of target individual
    filter(id == targ_id) %>%
    arrange(event_date)
  
  for(i in 1:nrow(ind_data)){ # assign herds to tibble column
    herd_id$herd[herd_id$year == ind_data$next_bio_year[i]] <- ind_data$population[i]
  }
  
  first_obs <- which.min(is.na(herd_id$herd))
  if(first_obs > 1){ # fill back one year, didn't come from nowhere
    herd_id$herd[first_obs - 1] <- herd_id$herd[first_obs]
  }
  
  for(i in first_obs:nrow(herd_id)){
    if(is.na(herd_id$herd[i])){ # fill forward to end 
      herd_id$herd[i] <- herd_id$herd[i-1]
    }
  }
  out <- c(targ_id, herd_id$herd) # build vector output
  names(out) <- c("id", herd_id$year)
  return(out)
}

herd_tib <- bind_rows(map( # apply over all ids in full data
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
      ) # only want events where the individual had to be alive
    ) %>%
  group_by(id) %>%
  summarize(
    first_capture = min(event_date),
    last_capture = max(event_date) # pull first/last live sighting of each id
  ) %>%
  mutate(first_alive = case_when( # Assign to capture years (alive on Nov 1)
    month(first_capture) <= cap_window[2] ~ year(first_capture) - 1,
    T ~ year(first_capture) # this mutate assumes calves not observed before Nov
  )) %>%
  mutate(last_alive = case_when( # assign last known alive year
    month(last_capture) >= cap_window[1] ~ year(last_capture),
    T ~ year(last_capture) - 1
  )) %>%
  ungroup()

alive_tib <- matrix(NA, nrow = nrow(herd_tib), ncol = ncol(herd_tib))

for(i in 1:nrow(herd_tib)){
  
  alive_tib[i,1] <- herd_tib$id[i] # share ids between tibs
  
  if(herd_tib$id[i] %in% f_l_data$id){ # do we have i in herd data?
    # assign first year, adjust if before study start
    first_year <- f_l_data[f_l_data$id == herd_tib$id[i],]$first_alive
    if(first_year < start_year){first_year <- start_year} 
    # assign last year, adjust if before study start
    last_year  <- f_l_data[f_l_data$id == herd_tib$id[i],]$last_alive
    if(last_year < start_year){last_year <- start_year} 
    # fill in ones for all years between first and last year
    # +2 because 1-1=0 and id takes a column
    alive_tib[i,(first_year:last_year) - start_year + 2] <- 1
  }
}

# get the tibble version together and sorted
alive_tib <- alive_tib %>% as_tibble()
names(alive_tib) <- c("id", as.character(start_year:(end_year+1)))
alive_tib <- alive_tib %>%
  arrange(id) %>%
  mutate_at(2:ncol(.), as.numeric)

# use known deaths to fill in trailing zeros
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
    alive_tib[i,death_column:ncol(alive_tib)] <- 0
  }
}

# use known births to fill in leading zeros and ones
birth_dates <- full_data %>%
  group_by(id) %>%
  summarise(birth_date = min(birth_date)) %>%
  ungroup() %>%
  arrange(id) %>%
  mutate(birth_year = case_when(
    month(birth_date) >= cap_window[1] ~ year(birth_date) + 1,
    T ~ year(birth_date)
  ))

birth_column <- rep(NA, nrow(alive_tib))
first_column <- rep(NA, nrow(alive_tib))

# This commented stuff doesn't work because of numeric issues. Grumble.
alive_tib <- alive_tib %>%
  select(-id) %>%
  as.matrix()
for(i in 1:nrow(alive_tib)){
  if(sum(alive_tib[i,], na.rm = T) > 0){
    first_column[i] <- min(which(alive_tib[i,] == 1))
    if(!is.na(birth_dates$birth_year[i])){
      # assertthat::assert_that(birth_dates$id[i] == alive_tib$id[i])
      birth_column[i] <- birth_dates$birth_year[i] - start_year + 1
      if(birth_column[i] < 1){birth_column[i] <- 1}
      if((first_column[i] - birth_column[i]) > 0){
        alive_tib[i,birth_column[i]:first_column[i]] <- 1
        # all of this is hacked together to sidestep errors, best of luck
        # untangling this mess future Kenneth
        # Sincerely,
        # That jerk from the past
      }
      alive_tib[i,(1:(birth_column[i] - 1))] <- 0
    }
  }
}

alive_tib <- alive_tib %>%
  as_tibble() %>%
  mutate(id = birth_dates$id) %>%
  select(id, 1:(ncol(.)-1))
  # mutate_at(vars(-id), funs(as.numeric))


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
    if(any(tmp$Class == "CA" | tmp$Class == "NC", na.rm = T)){ # if its captured as a calf
      birth_year <- year(min(tmp$event_date)) - (month(min(tmp$event_date)) < 5)
      out_date <- ymd(paste0(birth_year, "-06-01")) # use most recent June 1st
    }else{ # If there is no age information
      out_date <- NA # leave it blank
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
# okay. But calf session basically birth year in occasion units.
# age_data <- age_data %>%
#   mutate(calf_session = case_when(
#     calf_session > 0   ~ calf_session,
#     T ~ 0
#   ))

age_tib <- matrix(NA, nrow = nrow(age_data), ncol = ncol(herd_tib) - 1)
# for(i in 1:nrow(age_data)){
#   age_tib[i, age_data$calf_session[i]] <- 1
# }
# for(i in 1:nrow(age_tib)){
#   if(all(is.na(age_tib[i,]))){
#     age_tib[i,] <- 0
#   }else{
#     age_tib[i, (1 + which.min(age_tib[i,] == 1)):ncol(age_tib)] <- 0
#   }
# }
c_occ <- age_data$calf_session
for(i in 1:nrow(age_data)){
  if(is.na(c_occ[i])){
    age_tib[i,] <- NA
  }else{
    if(c_occ[i] < 1){
      age_tib[i,] <- 2:(ncol(age_tib) + 1)
    }else{
      age_tib[i,c_occ[i]:ncol(age_tib)] <- 1
      age_tib[i,c_occ[i]:ncol(age_tib)] <- cumsum(age_tib[i,c_occ[i]:ncol(age_tib)])
      age_tib[i,1:(c_occ[i]-1)] <- 0
    }
  }
}

age_tib <- age_data %>%
  select(id) %>%
  bind_cols(as_tibble(age_tib))
names(age_tib) <- c("id", as.character(start_year:(end_year + 1)))

#Capture_history================================================================
# You've checked everything above this as of 10:39 on 12 Sep 2022
cap_data <- full_data %>%
  # animals we know are alive
  filter(CaptureMethod %in% c(
    "TRAPPED", "HANDLD", "NETTED", "DARTED", "CLOVER", "MOVEDA",  "OBSRVD")) %>%
  # events that happen after the survey start date
  filter(event_date >= ymd(paste(start_year, cap_window[1], "1"))) %>%
  mutate(cap_month = month(event_date)) %>%
  # this is subsetting to specific months
  # we want that for capture histories, but not known alive histories
  filter(cap_month >= cap_window[1] | cap_month <= cap_window[2]) %>%
  mutate(session = case_when(
    cap_month >= cap_window[1] ~ year(event_date) - start_year + 1,
    cap_month <= cap_window[2] ~ year(event_date) - start_year
  )) %>%
  arrange(id, session)

cap_tib <- matrix(
  0, 
  nrow = nrow(herd_tib), 
  ncol = ncol(alive_tib) - 1) %>%
  as_tibble() %>%
  bind_cols(herd_tib$id, .) 
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
  explainer = "cap is capture history for js -- 
  age is year that id was a calf -- 
  sex is sex of id -- 
  liv is 1 for individual known alive, 0 for dead, NA for unknown -- 
  hrd is which herd id belonged to on 11-1 of each year -- 
  hnt is a 1 for the year an elk was removed from the population (hunting) --"
)

save(elk_data, file = elk_data_location)

rm(list = ls())
