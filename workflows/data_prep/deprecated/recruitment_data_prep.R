# Manually pull data from capture/handling summaries:
# This csv has the incoming counts from main from the capture/handling summaries
summary_count_data <- read_csv("data//min_n_handle_summaries.csv") %>%
  mutate(yr = year) %>%
  select(yr, fe_ad, ca) %>%
  filter(yr != "1991") %>%
  filter(yr != "1992") %>%
  filter(yr != "1993") %>%
  filter(yr != "1994") %>%
  filter(yr != "1995") %>%
  filter(yr != "1997")

dtlist <- list(
  data = list(
    n_calf = summary_count_data$ca,
    n_cow  = summary_count_data$fe_ad
  ),
  constants = list(
    n_years = nrow(summary_count_data),
    years = summary_count_data$yr
  ),
  initial_values = list(R = rep(0.5, nrow(summary_count_data)))
)

saveRDS(dtlist, "data//recruitment_data.rds")  
rm(list = ls())
