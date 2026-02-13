require(rjags); require(tidyverse)
rs_raw <- readRDS("results//fbipm_rslt_06feb2026_full_mean.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows()

locolor <- "#00a"
ntcolor <- "#a40"
ntdt <- rs_raw %>%
  select(grep("Ntot", names(.))) 
names(ntdt) <- as.character(1988:2023)
ntplot <- ntdt %>%
  pivot_longer(1:ncol(.), names_to = "yr", values_to = "val") %>%
  group_by(yr) %>%
  summarise(
    lci = quantile(val, 0.025),
    mci = quantile(val, 0.5),
    uci = quantile(val, 0.975)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = as.numeric(yr), y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = ntcolor, alpha = 0.2) +
  geom_line(color = ntcolor) +
  theme_classic() +
  xlab(element_blank()) + ylab("N Individuals") + 
  labs(title = "a - Elk Abundance") + 
  theme(axis.text.x = element_blank()) +
  xlim(c(1988, 2023))

lamplot <- readRDS("results//full_lambda_tibble_best_mean.rds") %>%
  select(yrx, Lo) %>%
  group_by(yrx) %>%
  summarise(
    lci = quantile(Lo, 0.025),
    mci = quantile(Lo, 0.5),
    uci = quantile(Lo, 0.975)
  ) %>% 
  ungroup() %>%
  mutate(yr = as.numeric(yrx) + 1987) %>%
  ggplot(aes(x = as.numeric(yr), y = mci)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), fill = locolor, alpha = 0.2) +
  geom_line(color = locolor) +
  geom_hline(yintercept = 1, linetype = "11", color = "#999") +
  theme_classic() +
  xlab("Year") + ylab("Lambda") + labs(title = "b - Population Growth") +
  xlim(c(1988, 2023))

cowplot::plot_grid(
  ntplot, lamplot,
  align = "v",
  ncol = 1)
ggsave(
  "figures/n_lambda_plot.png",
  width = 8,
  height = 8,
  units = "cm",
  dpi = 600,
  scale = 1.5)
