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

misc_data <- readRDS("data//misc_data.rds")
counts <- misc_data$count %>%
  filter(sex == "f" & age == "ad") %>%
  select(yr, count)
yr_range <- 1988:2023
#Capture Probability============================================================
sd_p <- rslt %>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(yr = as.numeric(yr)) %>%
  filter(yr != 2023) %>%
  full_join(counts) %>%
  filter(!is.na(count)) %>%
  mutate(est = count / val) %>%
  group_by(yr) %>%
  summarize(est_sd = sd(est)) %>%
  mutate(method = "Posterior") %>%
  select(yr, method, est_sd)
  
sd_afcount <- readRDS("s2//results//ipmrs_27sep2025_null.rds") %>%
  map(as_tibble) %>%
  bind_rows() %>%
  pull(sd_afcount)
sd_afcount %>% quantile(c(0.025, 0.5, 0.975))

sd_nb <- rslt %>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(yr_range + 1)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "val") %>%
  filter(yr %in% as.character(yr_range[-length(yr_range)] + 1)) %>%
  mutate(yr = as.numeric(yr)) %>%
  filter(yr != 2023) %>%
  group_by(yr) %>%
  summarize(pm = quantile(val, 0.5)) %>%
  full_join(counts) %>%
  filter(!is.na(count)) %>%
  mutate(est_sd = sqrt(count*(1-pm)/(pm^2))) %>%
  mutate(method = "NB") %>%
  select(yr, method, est_sd)

bind_rows(sd_nb, sd_p) %>%
  rm_bad_yrs() %>%
  ggplot(aes(x = yr, y = est_sd, color = method, shape = method)) +
  geom_line() + geom_point() +
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab("Year") + ylab("Estimated SD")
ggsave("figures//count_sd_method_comp.png", dpi = 600, units = "cm",
       height = 7, width = 8.5)
