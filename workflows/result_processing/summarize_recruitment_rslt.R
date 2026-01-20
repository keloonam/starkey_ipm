yrs <- readRDS("data//recruitment_data.rds")$constants$years
logit <- function(x){log(x/(1-x))}
rslt_summ <- readRDS("results//recruitment_rslt.rds") %>%
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


saveRDS(rslt_summ, file = "results//recr_rslt_summ.rds")
