require(tidyverse)
edt <- readRDS("results//lambda_elasticity.rds")
sdt <- readRDS("results//lambda_sensitivity.rds")
adt <- readRDS("results//annual_L_effects.rds")
edt$wvcv %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975)
  ) %>%
  ggplot(aes(x = param, y = mci)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = 0.5))

sdt %>%
  # filter(param %in% c("cg", "dd", "wt", "cm", "dm", "wm")) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(val, 0.025),
    mci = quantile(val, 0.5),
    uci = quantile(val, 0.975)
  ) %>%
  ggplot(aes(x = param, y = mci)) +
  geom_pointrange(aes(ymin = lci, ymax = uci))

ann_dt <- apply(adt, 2:3, quantile, probs = 0.5) %>%
  as_tibble() %>%
  mutate(
    Puma = cg + cm, 
    Elk = dd + dm,
    Ratio = nfp + ncp,
    Recruitment = ppyr + snyr,
    Survival = scyr + sfyr,
    noise = ppyr + snyr + scyr + sfyr,
    Climate = wm + wt,
    Harvest = sh1 + sh2) %>%
  mutate(yr = 2:(nrow(.)+1)) %>%
  pivot_longer(1:(ncol(.)-1), names_to = "param", values_to = "val")
ann_dt %>%
  filter(param %in% c(
    "ppyr", "sfyr", "scyr", "snyr"
    )) %>%
  ggplot(aes(x = yr, y = val, fill = param)) +
  geom_col()

ann_dt %>%
  filter(param %in% c(
    "ppyr", "sfyr", "scyr", "snyr"
  )) %>%
  group_by(param) %>%
  summarise(mn = mean(abs(val)))
