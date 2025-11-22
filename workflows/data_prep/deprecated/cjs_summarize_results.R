# Process CJS results down to summary stats on the logit scale
# This format is needed for the IPM, but these objects are not quite IPM ready

#Functions======================================================================
logit <- function(x){
  log(x/(1-x)) %>% return()
}

rm_bad_yrs <- function(x){
  x %>%
    filter(yr != "2017") %>%
    filter(yr != "2020") %>%
    return()
}

#Load the results===============================================================
rslt <- readRDS("results//cjs_rslt.rds") %>%
  map(as_tibble) %>%
  bind_rows()

#Capture Probability============================================================
fp_logit <- rslt %>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bpf")
mp_logit <- rslt %>%
  select(grep("prob_am", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bpm")

#Survival=======================================================================
fs_logit <- rslt %>%
  select(grep("survival_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsf")
ms_logit <- rslt %>%
  select(grep("survival_am", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsm")
cs_logit <- rslt %>%
  select(grep("survival_ca", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(val = logit(val)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val)) %>%
  ungroup() %>%
  mutate(var = "bsc")

#Cleanup Steps==================================================================
rslt_summary <- list(
  sf_logit = fs_logit,
  sm_logit = ms_logit,
  sc_logit = cs_logit,
  pf_logit = fp_logit,
  pm_logit = mp_logit
)


rslt_summary <- map(rslt_summary, rm_bad_yrs)


saveRDS(rslt_summary, "results//cjs_summary.rds")
rm(list = ls()[-which(ls() == "yr_range")])
