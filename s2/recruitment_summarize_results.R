# Prepare summary stats of recruitment results on the log scale
# This is the format needed for the IPM but these objects are not IPM ready yet

#Functions======================================================================
logit <- function(x){
  log(x/(1-x)) %>% return()
}

#Load Results===================================================================
yrs <- readRDS("data//recruitment_data.rds")$constants$years
rs <- readRDS("results//recruitment_rslt_10sep2024.rds") %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("R", names(.))) %>%
  set_names(as.character(yrs)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "r") %>%
  group_by(yr) %>%
  summarise(mn = mean(r), tau = 1/sd(r)^2) %>%
  ungroup() %>%
  mutate(yr = as.numeric(yr))

#Cleanup========================================================================
saveRDS(rs, "s2//recruitment_summary.rds")
rm(list = ls())
