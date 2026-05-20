require(tidyverse)
save_file_name <-  "figures//elasticity_comparison_SI.png"
tag <- "best_lgst"
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
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    prc, pr2, pry, pro, prp,
    sh1, sh2,
    ppyr, snyr, scyr, sfyr
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95)
  ) %>%
  ungroup() %>%
  mutate(prm = case_when(
    param == "dd" ~ "l dd",
    param == "wt" ~ "m wt",
    param == "wm" ~ "n wm",
    param == "cg" ~ "o cg",
    param == "ppyr" ~ "a pregnancy",
    param == "snyr" ~ "b neonate survival",
    param == "scyr" ~ "c calf survival",
    param == "sfyr" ~ "d adult survival",
    param == "sh1" ~ "e calf harvest",
    param == "sh2" ~ "f adult harvest",
    param == "prc" ~ "g calf",
    param == "pr2" ~ "h yearling",
    param == "pry" ~ "i young",
    param == "prp" ~ "j prime",
    param == "pro" ~ "k old"
  )) %>%
  mutate(tag = tag)

tag <- "best_mean"
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
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    prc, pr2, pry, pro, prp,
    sh1, sh2,
    ppyr, snyr, scyr, sfyr
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95)
  ) %>%
  ungroup() %>%
  mutate(prm = case_when(
    param == "dd" ~ "l dd",
    param == "wt" ~ "m wt",
    param == "wm" ~ "n wm",
    param == "cg" ~ "o cg",
    param == "ppyr" ~ "a pregnancy",
    param == "snyr" ~ "b neonate survival",
    param == "scyr" ~ "c calf survival",
    param == "sfyr" ~ "d adult survival",
    param == "sh1" ~ "e calf harvest",
    param == "sh2" ~ "f adult harvest",
    param == "prc" ~ "g calf",
    param == "pr2" ~ "h yearling",
    param == "pry" ~ "i young",
    param == "prp" ~ "j prime",
    param == "pro" ~ "k old"
  )) %>%
  mutate(tag = tag)

tag <- "best_norm"
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
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    prc, pr2, pry, pro, prp,
    sh1, sh2,
    ppyr, snyr, scyr, sfyr
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95)
  ) %>%
  ungroup() %>%
  mutate(prm = case_when(
    param == "dd" ~ "l dd",
    param == "wt" ~ "m wt",
    param == "wm" ~ "n wm",
    param == "cg" ~ "o cg",
    param == "ppyr" ~ "a pregnancy",
    param == "snyr" ~ "b neonate survival",
    param == "scyr" ~ "c calf survival",
    param == "sfyr" ~ "d adult survival",
    param == "sh1" ~ "e calf harvest",
    param == "sh2" ~ "f adult harvest",
    param == "prc" ~ "g calf",
    param == "pr2" ~ "h yearling",
    param == "pry" ~ "i young",
    param == "prp" ~ "j prime",
    param == "pro" ~ "k old"
  )) %>%
  mutate(tag = tag)

tag <- "best_odfw"
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
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    prc, pr2, pry, pro, prp,
    sh1, sh2,
    ppyr, snyr, scyr, sfyr
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95)
  ) %>%
  ungroup() %>%
  mutate(prm = case_when(
    param == "dd" ~ "l dd",
    param == "wt" ~ "m wt",
    param == "wm" ~ "n wm",
    param == "cg" ~ "o cg",
    param == "ppyr" ~ "a pregnancy",
    param == "snyr" ~ "b neonate survival",
    param == "scyr" ~ "c calf survival",
    param == "sfyr" ~ "d adult survival",
    param == "sh1" ~ "e calf harvest",
    param == "sh2" ~ "f adult harvest",
    param == "prc" ~ "g calf",
    param == "pr2" ~ "h yearling",
    param == "pry" ~ "i young",
    param == "prp" ~ "j prime",
    param == "pro" ~ "k old"
  )) %>%
  mutate(tag = tag)

tag <- "best_rcns"
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
  select(
    stpn, tot, 
    cg, dd, wm, wt, 
    prc, pr2, pry, pro, prp,
    sh1, sh2,
    ppyr, snyr, scyr, sfyr
  ) %>%
  pivot_longer(3:ncol(.), names_to = "param", values_to = "val") %>%
  mutate(value = abs(val) / tot) %>%
  group_by(param) %>%
  summarise(
    lci = quantile(value, 0.05),
    mci = quantile(value, 0.5),
    uci = quantile(value, 0.95)
  ) %>%
  ungroup() %>%
  mutate(prm = case_when(
    param == "dd" ~ "l dd",
    param == "wt" ~ "m wt",
    param == "wm" ~ "n wm",
    param == "cg" ~ "o cg",
    param == "ppyr" ~ "a pregnancy",
    param == "snyr" ~ "b neonate survival",
    param == "scyr" ~ "c calf survival",
    param == "sfyr" ~ "d adult survival",
    param == "sh1" ~ "e calf harvest",
    param == "sh2" ~ "f adult harvest",
    param == "prc" ~ "g calf",
    param == "pr2" ~ "h yearling",
    param == "pry" ~ "i young",
    param == "prp" ~ "j prime",
    param == "pro" ~ "k old"
  )) %>%
  mutate(tag = tag)




bind_rows(lgst_dt, mean_dt, norm_dt, odfw_dt, rcns_dt) %>%
  ggplot(aes(x = prm, y = mci, color = tag, shape = tag)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge(width = .5), size = 0.5) +
  theme_classic() +
  scale_x_discrete(
    labels = c(
      `l dd` = "Density",
      `m wt` = bquote("SPEI"[t]),
      `n wm` = bquote("SPEI"[t-1]),
      `o cg` = "Pumas",
      `a pregnancy` = "Pregnancy",
      `b neonate survival` = "S (pre-wean)",
      `c calf survival` = "S (post-wean)",
      `d adult survival` = "S (adult)",
      `e calf harvest` = "H (calf)",
      `f adult harvest` = "H (adult)",
      `g calf` = "PN (calf)",
      `h yearling` = "PN (yearling)",
      `i young` = "PN (young)",
      `j prime` = "PN (prime)",
      `k old` = "PN (old)"
    ),
    guide = guide_axis(n.dodge = 2)
  ) +
  scale_color_discrete(
    labels = c(
      best_lgst = "Logistic Growth",
      best_mean = "Mean Index",
      best_norm = "Random Normal",
      best_odfw = "ODFW Estimate",
      best_rcns = "Reconstruction"
    ),
    palette = c(
      "#0072B2", 
      "#D55E00",
      "#E69F00",
      "#009E73",
      "#800080"
    )) +
  scale_shape_discrete(
    labels = c(
      best_lgst = "Logistic Growth",
      best_mean = "Mean Index",
      best_norm = "Random Normal",
      best_odfw = "ODFW Estimate",
      best_rcns = "Reconstruction"
    )) +
  xlab(NULL) + ylab("Contribution") +
  labs(color = element_blank(), shape = element_blank()) +
  theme(legend.position = "bottom") 
ggsave(
  filename = paste0(save_file_name),
  width = 18,
  height = 10,
  units = "cm",
  dpi = 600
)
