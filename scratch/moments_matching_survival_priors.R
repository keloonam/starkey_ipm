require(dplyr); require(purrr); require(rjags); require(tidyr)
logit <- function(p){
  out <- log(p/(1-p)) 
  return(out)
}
ilogit <- function(x){
  out <- exp(x)/(1+exp(x)) 
  return(out)
}


# female 
rlogis(1000000, 2.528456, 0.3831403) |> ilogit() |> hist(breaks = 100)
# male
rlogis(1000000, 1.466578, 0.3419404) |> ilogit() |> hist(breaks = 100)

1/sqrt(3)
3^(-1/2)


(logit(0.92) - logit(0.95)) / 1.95 * pi * 3^(-1/2)

logit(0.95)
ilogit(2)

fs_real <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_af")) %>%
  pivot_longer(cols = everything(), names_to = "yr", values_to = "val") %>%
  mutate(nv = logit(val)) %>%
  group_by(yr) %>%
  summarise(myr = mean(nv), syr = sd(nv)) %>%
  pull(myr)

ms_real <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_am")) %>%
  pivot_longer(cols = everything(), names_to = "yr", values_to = "val") %>%
  mutate(nv = logit(val)) %>%
  group_by(yr) %>%
  summarise(myr = mean(nv), syr = sd(nv)) %>%
  pull(myr)

mean(ms_real); sd(ms_real); sd_to_scale(sd(ms_real))
mean(fs_real); sd(fs_real); sd_to_scale(sd(fs_real))

sd_to_scale <- function(s){
  s * sqrt(3) / pi
}
sd_to_precision <- function(s){
  1/(s^2)
}

fs_real_2 <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_af")) %>%
  as.matrix() %>%
  logit()
apply(fs_real_2, 1, sd) %>% mean(); apply(fs_real_2, 1, sd) %>% sd() %>% sd_to_precision()


ms_real_2 <- rslt %>%
  map(as_tibble) %>%
  bind_rows() %>%
  select(contains("survival_am")) %>%
  as.matrix() %>%
  logit()
apply(ms_real_2, 1, sd) %>% mean(); apply(ms_real_2, 1, sd) %>% sd() %>% sd_to_precision()


rnorm(1000000, 0.83, sqrt(1/45.64)) %>% hist()
