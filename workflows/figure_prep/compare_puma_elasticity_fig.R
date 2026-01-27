require(tidyverse)
tag <- "bst_lgst"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
lgst_dt <- edt$wvcv %>%
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
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  mutate(tag = tag)

tag <- "bst_mean"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
mean_dt <- edt$wvcv %>%
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
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  mutate(tag = tag)

tag <- "bst_norm"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
norm_dt <- edt$wvcv %>%
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
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  mutate(tag = tag)

tag <- "bst_odfw"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
odfw_dt <- edt$wvcv %>%
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
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  mutate(tag = tag)

tag <- "bst_rcns"
edt_name <- paste0("results//lambda_elasticity_", tag, ".rds")
sdt_name <- paste0("results//lambda_sensitivity_", tag, ".rds")
edt <- readRDS(edt_name)
sdt <- readRDS(sdt_name)
total_variation_by_step <- edt$wvcv %>%
  mutate(val = abs(value)) %>%
  summarise(tot = sum(val), .by = stpn)
rcns_dt <- edt$wvcv %>%
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
  mutate(prm = case_when(
    param == "dd" ~ "d dd",
    param == "wt" ~ "e wt",
    param == "wm" ~ "f wm",
    param == "cg" ~ "g cg",
    param == "annual_variation" ~ "b demographic variation",
    param == "harvest" ~ "a harvest",
    param == "age_structure" ~ "c age distribution"
  )) %>%
  mutate(tag = tag)




bind_rows(lgst_dt, mean_dt, norm_dt, odfw_dt, rcns_dt) %>%
  ggplot(aes(x = prm, y = mci, color = tag)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = .4)) +
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
  scale_color_discrete(
    labels = c(
      bst_lgst = "Logistic Growth",
      bst_mean = "Index Mean",
      bst_norm = "Normal Draw",
      bst_odfw = "ODFW Estimate",
      bst_rcns = "Reconstruction"
  )) +
  xlab(NULL) + ylab("Contribution") +
  labs(color = "Puma Covariate")
ggsave(
  filename = paste0("figures//lambda_cont_puma_comparison.png"),
  width = 18,
  height = 10,
  units = "cm",
  dpi = 600
)
