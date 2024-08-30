build_herd_history <- function(targ_id, full_data, start_year, end_year){
  # Assigns an individual to a heard on each Nov 1 after its first observation
  
  herd_id <- tibble( # empty single-individual tibble
    year = start_year:end_year,
    herd = NA
  )
  
  ind_data <- full_data %>% # get events of target individual
    filter(id == targ_id) %>%
    arrange(event_dt)
  
  for(i in 1:nrow(ind_data)){ # assign herds to tibble column
    # This will overwrite the herd each time a year repeats
    # Because it is pre-sorted, the most recent herd is recorded
    herd_id$herd[herd_id$year == ind_data$next_year[i]] <- ind_data$pop[i]
  }
  
  first_obs <- which.min(is.na(herd_id$herd)) # records position of first non-NA
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

build_known_alive_tibble <- function(targ_id, full_data, start_yr, end_yr){
  
  # get events of target individual
  ind_data <- full_data %>% 
    filter(id == targ_id) %>%
    arrange(event_dt)
  
  # pull relevant dates
  birth_date <- min(ind_data$birth_dt)
  death_date <- max(ind_data$death_dt)
  event_frst <- min(ind_data$event_dt)
  event_last <- max(ind_data$event_dt)
  
  # convert to years then indices, rather than dates
  birth_t <- year(birth_date) + (month(birth_date) > 10) - start_yr + 1
  death_t <- year(death_date) + (month(death_date) > 10) - start_yr + 1
  first_t <- year(event_frst) - (month(event_frst) < 11) - start_yr + 1
  final_t <- year(event_last) - (month(event_last) < 11) - start_yr + 1
  
  life_history <- rep(NA, end_yr - start_yr + 1) # Prepare empty vector
  life_history[first_t:final_t] <- 1 # alive between captures
  if(!is.na(birth_t) & birth_t >= 1){
    life_history[0:(birth_t - 1)] <- 0 # "dead" before birth
    life_history[birth_t:final_t] <- 1 # alive between birth and last capture
  }
  if(!is.na(death_t)){
    life_history[first_t:(death_t - 1)] <- 1 # alive between first cap and death
    life_history[death_t:length(life_history)] <- 0 # dead after death
  }
  
  out <- c(targ_id, life_history) # combine with id
  names(out) <- c("id", start_yr:end_yr) # name the "columns"
  return(out)
}

build_age_history <- function(targ_id, full_data, start_yr, end_yr){
  ind_data <- full_data %>%  # subset down to data for one individual
    filter(id == targ_id) %>%
    arrange(event_dt) %>%
    select(id, Class, event_dt, birth_dt, death_dt)
  
  suppressWarnings(
    birth_t <- ind_data %>% # pull the birth date
      mutate(birth_t = year(birth_dt) + (month(birth_dt)>10) - start_yr + 1) %>%
      pull(birth_t) %>% min(na.rm = T) 
  )
  suppressWarnings(
    if(is.na(birth_t)){ # if no birth date
      birth_t <- ind_data %>% # infer from capture when it was a calf
        filter(Class %in% c("CA", "NC")) %>%
        mutate(birth_t = year(event_dt) - (month(event_dt)<5)- start_yr + 1) %>%
        pull(birth_t) %>% min()
    }
  )

  age_history <- rep(NA, length(start_yr:end_yr)) # empty vector
  if(birth_t != Inf){
    if(birth_t < 1){ # if they were born before the first year
      age_history[1] <- abs(birth_t) + 1 # calculate their age in year 1
      age_history[2:length(age_history)] <- 1
    }else{
      age_history[0:(birth_t - 1)] <- 0 # add 0s before birth
      age_history[birth_t:length(age_history)] <- 1 # 1s from birth onward
    }
    age_history <- cumsum(age_history) # age is cumulative sum
    age_history[age_history == 0] <- NA # switch pre-birth back to NA
  }
  
  out <- c(targ_id, age_history) # combine with id
  names(out) <- c("id", start_yr:end_yr) # name the "columns"
  return(out)
}

build_capture_history <- function(full_data, start_yr, end_yr, cap_window){
  cap_data <- full_data %>%
    # animals we know are alive
    filter(CaptureMethod %in% c(
      "TRAPPED", "HANDLD", "NETTED", "DARTED", "CLOVER", "MOVEDA", "OBSRVD"
    )) %>%
    # events that happen after the survey start date
    filter(event_dt >= ymd(paste(start_yr, cap_window[1], "1"))) %>%
    mutate(cap_month = month(event_dt)) %>%
    # this is sub-setting to specific months
    # we want that for capture histories, but not known alive histories
    filter(cap_month >= cap_window[1] | cap_month <= cap_window[2]) %>%
    mutate(session = case_when(
      cap_month >= cap_window[1] ~ year(event_dt) - start_yr + 1,
      cap_month <= cap_window[2] ~ year(event_dt) - start_yr
    )) %>%
    arrange(id, session)
  
  id_col <- tibble(id = unique(full_data$id))
  
  cap_tib <- matrix(
    0, 
    nrow = length(unique(full_data$id)), 
    ncol = length(start_yr:end_yr),
    dimnames = list(
      NULL, c(as.character(start_yr:(end_yr)))
    )) %>%
    as_tibble() %>%
    bind_cols(id_col, .) 
  
  for(i in 1:nrow(cap_data)){
    cap_tib[cap_tib$id == cap_data$id[i], cap_data$session[i] + 1] <- 1
  }
  return(cap_tib)
}

build_harvest_history <- function(full_data, start_yr, end_yr){
  hunt_data <- full_data %>%
    filter(CaptureMethod == "HUNTED") %>%
    mutate(cap_month = month(event_dt)) %>%
    mutate(session = case_when(
      cap_month >= 11 ~ year(event_dt) - start_yr + 2,
      T ~ year(event_dt) - start_yr + 1
    )) %>%
    arrange(id, session)
  
  id_col <- tibble(id = unique(full_data$id))
  
  hunt_tib <- matrix(
    0, 
    nrow = length(unique(full_data$id)), 
    ncol = length(start_yr:end_yr),
    dimnames = list(
      NULL, c(as.character(start_yr:(end_yr)))
    )) %>%
    as_tibble() %>%
    bind_cols(id_col, .) 
  
  for(i in 1:nrow(hunt_data)){
    hunt_tib[hunt_tib$id == hunt_data$id[i], hunt_data$session[i] + 1] <- 1
  }
  return(hunt_tib)
}
