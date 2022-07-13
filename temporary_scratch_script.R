cap_data <- full_data %>%
  # animals we know are alive
  filter(CaptureMethod %in% c(
    "TRAPPED", "HANDLD", "NETTED", "DARTED", "CLOVER", "MOVEDA",  "OBSRVD",
    "HUNTED", "MRT-TRAP")) %>%
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
  nrow = length(unique(cap_data$id)), 
  ncol = ncol(alive_tib) - 1) %>%
  as_tibble() %>%
  bind_cols(unique(cap_data$id), .) 
names(cap_tib) <- c("id", as.character(start_year:(end_year + 1)))

for(i in 1:nrow(cap_data)){
  cap_tib[cap_tib$id == cap_data$id[i], cap_data$session[i] + 1] <- 1
}