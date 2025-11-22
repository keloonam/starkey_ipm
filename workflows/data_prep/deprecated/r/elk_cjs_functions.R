# Functions supporting elk cjs data cleaning
# Kenneth Loonam
# March 2021

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
        }else{
          out_date <- NA
        }
      }
    }
  }
  
  out <- tibble(id = id_x, birth_date = out_date)
}
