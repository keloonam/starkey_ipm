require(tidyverse); require(rjags)

mcmc_to_tibble <- function(raw_file, session, bayes){
  raw_file %>%
    map(as_tibble) %>%
    bind_rows() %>%
    mutate(session = session) %>%
    mutate(bayes = bayes) %>%
    return()
}

pull_demographic <- function(x, ch_pattern){
  x %>%
    select(c(session, bayes, grep(ch_pattern, names(.)))) %>%
    setNames(c("session", "bayes", as.character(1:ncol(.)))) %>%
    pivot_longer(cols = 3:ncol(.), names_to = "yr_index", values_to = "val") %>%
    group_by(session, bayes, yr_index) %>%
    summarize(
      lci = quantile(val, 0.025),
      mci = quantile(val, 0.5),
      uci = quantile(val, 0.975),
      sd = sd(val)
    ) %>%
    ungroup() %>%
    mutate(yr_index = as.numeric(yr_index)) %>%
    return()
}

file_names <- paste0("fb_tests//results//", list.files("fb_tests//results"))
argument_list <- list(
  raw_file = map(file_names, readRDS),
  session = c(1:6, 1:6),
  bayes = c(rep("fb", 6), rep("pb", 6))
)

frs <- pmap(.l = argument_list, .f = mcmc_to_tibble) %>%
  bind_rows()

recruitment <- frs %>%
  pull_demographic("R\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd)

ggplot(recruitment, aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymax = uci, ymin = lci), position = position_dodge2(0.5))
  
af_survival <- frs %>%
  pull_demographic("SF\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd)
ggplot(af_survival, aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymax = uci, ymin = lci), position = position_dodge2(0.5))


recruitment %>% 
  filter(yr < 2018) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  pivot_wider(names_from = bayes, values_from = c(lci, mci, uci, sd)) %>%
  ggplot(aes(x = mci_pb, y = mci_fb)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

ca_survival <- frs %>%
  pull_demographic("SC\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd)
ggplot(ca_survival, aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymax = uci, ymin = lci), position = position_dodge2(0.5))

ntot <- frs %>%
  pull_demographic("NC\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd)
ggplot(ntot, aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymax = uci, ymin = lci), position = position_dodge2(0.5))
