logit <- function(x){log(x/(1-x)) %>% return()}
ilogit <- function(x){
  (1/(1+exp(-x))) %>% return()
}
rdt <- new_res %>% map(as_tibble) %>%
  bind_rows() %>%
  select(grep("R\\[", names(.))) %>%
  set_names(as.character(1989:2023)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "R")
rdt <- rdt %>%
  group_by(yr) %>%
  summarise(R = mean(R)) %>%
  ungroup() %>%
  mutate(yr = as.numeric(yr))
load("results/ipm_result_21apr2023_R&S_pdi.Rdata")
old_rdt <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(grep("R\\[", names(.))) %>%
  set_names(as.character(1989:2021)) %>%
  pivot_longer(cols = 1:ncol(.), names_to = "yr", values_to = "R") %>%
  group_by(yr) %>%
  summarise(R = mean(R)) %>%
  ungroup() %>%
  mutate(yr = as.numeric(yr))
plot(old_rdt)
lines(rdt)
load("data//elk_ipm_data_07jun2024.rdata")
ipm_data$r_ratio[,4]
jags_data$est_br
(ilogit(jags_data$est_br[,3])[-c(26,27)] - ipm_data$r_ratio[,4]) %>% hist()
