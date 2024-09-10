logit <- function(x){
  log(x/(1-x)) %>% return()
}

rm_bad_yrs <- function(x){
  x %>%
    filter(yr != "2017") %>%
    filter(yr != "2020") %>%
    return()
}

rslt <- readRDS(results_file) %>%
  map(as_tibble) %>%
  bind_rows()

fp_logit <- rslt %>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(1989:2022)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
mp_logit <- rslt %>%
  select(grep("prob_am", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(1989:2022)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
fs_logit <- rslt %>%
  select(grep("survival_af", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(1989:2022)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
ms_logit <- rslt %>%
  select(grep("survival_am", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(1989:2022)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
cs_logit <- rslt %>%
  select(grep("survival_ca", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(1989:2022)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))

rslt_summary <- list(
  sf_logit = fs_logit,
  sm_logit = ms_logit,
  sc_logit = cs_logit,
  pf_logit = fp_logit,
  pm_logit = mp_logit
)

if(remove_bad_years == T){
  rslt_summary <- map(rslt_summary, rm_bad_yrs)
}

saveRDS(rslt_summary, "results//cjs_summary.rds")
rm(list = ls())