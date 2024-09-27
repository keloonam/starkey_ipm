require(dplyr); require(purrr)
logit <- function(x){
  log(x/(1-x)) %>% return()
}
ilogit <- function(x){
  (1/(1+exp(-x))) %>% return()
}
load("results//survival//cjs_rslt_13sep2022.Rdata")
ors <- rslt %>%
  map(as_tibble) %>%
  bind_rows()

combine_cjs_results <- function(x, y, v){
  out <- tibble(
    yr = c(1989:2016, 2018, 2019, 2021, 2022),
    mn = NA,
    tau = NA,
    var = v
  )
  for(i in 1:nrow(out)){
    if(out$yr[i] %in% c(1989:2016, 2018, 2019)){
      out$mn[i] <- x$mn[which(x$yr == out$yr[i])]
      out$tau[i] <- 1/(x$sd[which(x$yr == out$yr[i])]^2)
    }else{
      out$mn[i] <- y$mn[which(y$yr == out$yr[i])]
      out$tau[i] <- 1/(y$sd[which(y$yr == out$yr[i])]^2)
    }
  }
  return(out)
}

o_sf <- ors%>%
  select(grep("survival_af", names(.))) %>%
  set_names(as.character(1989:2022)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
o_sm <-ors%>%
  select(grep("survival_am", names(.))) %>%
  set_names(as.character(1989:2022)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
o_sc <- ors %>%
  select(grep("survival_ca", names(.))) %>%
  set_names(as.character(1989:2022)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
o_pf <- ors%>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(1989:2022)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
o_pm <- ors%>%
  select(grep("prob_am", names(.))) %>%
  set_names(as.character(1989:2022)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))

nrs <- readRDS("results//cjs_rslt_24sep2024.rds") %>%
  map(as_tibble) %>%
  bind_rows()
n_sf <- nrs%>%
  select(grep("survival_af", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
n_sm <-nrs%>%
  select(grep("survival_am", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
n_sc <- nrs%>%
  select(grep("survival_ca", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
n_pf <- nrs%>%
  select(grep("prob_af", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
n_pm <- nrs%>%
  select(grep("prob_am", names(.))) %>%
  set_names(as.character(1989:2024)) %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(mn = mean(val), sd = sd(val))
# with(n_sf, plot(ilogit(mn) ~ yr))
# with(o_sf, lines(ilogit(mn) ~ yr))

cleaned_rslt <- list(
  sf_logit = combine_cjs_results(o_sf, n_sf, "bsf"),
  sm_logit = combine_cjs_results(o_sm, n_sm, "bsm"),
  sc_logit = combine_cjs_results(o_sc, n_sc, "bsc"),
  pf_logit = combine_cjs_results(o_pf, n_pf, "bpf"),
  pm_logit = combine_cjs_results(o_pm, n_pm, "bpm")
)
saveRDS(cleaned_rslt, "s2//cjs_summary.rds")

rm(list = ls())
