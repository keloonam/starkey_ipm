# Prepare summary stats of recruitment results on the log scale
# This is the format needed for the IPM but these objects are not IPM ready yet

#Functions======================================================================
logit <- function(x){
  log(x/(1-x)) %>% return()
}

#Load Results===================================================================
yrs <- readRDS(recruitment_data)$constants$years
rs <- readRDS(results_file) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("R", names(.))) %>%
  set_names(as.character(yrs)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "r") %>%
  mutate(rb0 = logit(r)) %>%
  group_by(yr) %>%
  summarise(mn = mean(rb0), sd = sd(rb0)) %>%
  ungroup() %>%
  mutate(yr = as.numeric(yr))

#Cleanup========================================================================
saveRDS(rs, "results//recruitment_summary.rds")
rm(list = ls())
