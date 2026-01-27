require(tidyverse)
tag <- "bst_lgst"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
# adt <- readRDS("results//annual_L_effects.rds")
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
bar_graph_data <- edt$wvcv %>%
  pivot_wider(id_cols = stpn, names_from = param, values_from = value) %>%
  left_join(total_variation_by_step) %>%
  mutate(
    age_structure = prc + pr2 + pry + pro + prp,
    harvest = sh1 + sh2,
    annual_variation = ppyr + snyr + scyr + sfyr
  ) %>%
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    age_structure, harvest, annual_variation
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  ))

edt$wvcv %>%
  pivot_wider(id_cols = stpn, names_from = param, values_from = value) %>%
  left_join(total_variation_by_step) %>%
  mutate(
    age_structure = prc + pr2 + pry + pro + prp,
    harvest = sh1 + sh2,
    annual_variation = ppyr + snyr + scyr + sfyr
  ) %>%
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    age_structure, harvest, annual_variation
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.025),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.975)
  ) %>%
  ungroup() %>%
  # mutate(var_group = case_when(
  #   param %in% c("cg", "dd", "wm", "wt") ~ "a covariate",
  #   param == "age_structure" ~ "c age distribution",
  #   param =="harvest" ~ "d harvest",
  #   param == "annual_variation" ~ "b demographic variation"
  # )) %>%
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  ggplot(aes(x = prm, y = mci)) +
  geom_pointrange(aes(ymin = lci, ymax = uci)) +
  theme_classic() +
  scale_x_discrete(
    labels = c(
      `d dd` = "Density",
      `e wt` = bquote("SPEI"[t]),
      `f wm` = bquote("SPEI"[t-1]),
      `g cg` = "Pumas",
      `b demographic variation` = "Random Effects",
      `a harvest` = "Harvest",
      `c age distribution` = "Age Structure"
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  xlab(NULL) + ylab("Contribution")
ggsave(
  filename = paste0("figures//lambda_cont_", tag, ".png"),
  width = 12,
  height = 8,
  units = "cm",
  dpi = 600
)