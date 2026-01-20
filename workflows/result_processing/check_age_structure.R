require(tidyverse); require(rjags)

rslt <- readRDS("results//ipm_rslt_17dec2025_spei12.rds") %>%
  map(as.matrix) %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(chain = c(
    rep("a", nrow(.)/5),
    rep("b", nrow(.)/5),
    rep("c", nrow(.)/5),
    rep("d", nrow(.)/5),
    rep("e", nrow(.)/5)
    )) %>%
  mutate(mcmc_step = 1:nrow(.)) %>%
  select(grep("N\\[", names(.)))

ndt <- rslt %>%
  pivot_longer(1:ncol(.)) %>%
  filter(!grepl("SN", name)) %>%
  separate_wider_delim(cols = name, delim = "[", names = c("a","b")) %>%
  separate_wider_delim(cols = b, delim = "]", names = c("b", "c")) %>%
  separate_wider_delim(cols = b, delim = ",", names = c("age", "sex", "yr")) %>%
  mutate(age = as.numeric(age), sex = as.numeric(sex), yr = as.numeric(yr)) %>%
  mutate(N = value) %>%
  select(age, sex, yr, value)

pdt <- ndt %>%
  mutate(stpn = sort(rep(
    1:nrow(rslt), 
    length(unique(ndt$age)) * length(unique(ndt$yr)) * length(unique(ndt$sex))
    ))) %>%
  filter(sex == 1) %>%
  mutate(class = case_when(
    age == 1 ~ "calf",
    age <= 3 ~ "young",
    age < 14 ~ "prime",
    T ~ "old"
  )) %>%
  group_by(class, yr, stpn) %>%
  summarise(N = sum(value)) %>%
  ungroup() %>%
  pivot_wider(values_from = N, names_from = class) %>%
  mutate(ntot = calf + old + prime + young) %>%
  mutate(
    calf = calf / ntot,
    young = young / ntot,
    prime = prime / ntot,
    old = old / ntot
  ) %>%
  select(yr, calf, young, prime, old) %>%
  pivot_longer(2:ncol(.), names_to = "class", values_to = "N") %>%
  group_by(class, yr) %>%
  summarise(
    lci = quantile(N, 0.025),
    mci = quantile(N, 0.5),
    uci = quantile(N, 0.975)
  )

pdt %>%
  ggplot(aes(x = yr, y = mci, color = class)) +
  geom_line() +
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = class), alpha = 0.5)
