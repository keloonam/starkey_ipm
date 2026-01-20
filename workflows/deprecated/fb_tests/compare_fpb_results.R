#Environment====================================================================
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
      lci = quantile(val, 0.025, na.rm = T),
      mci = quantile(val, 0.5, na.rm = T),
      uci = quantile(val, 0.975, na.rm = T),
      sd = sd(val, na.rm = T)
    ) %>%
    ungroup() %>%
    mutate(yr_index = as.numeric(yr_index)) %>%
    return()
}

#Data===========================================================================

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
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024)
af_survival <- frs %>%
  pull_demographic("SF\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024)
am_survival <- frs %>%
  pull_demographic("SM\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024)
ca_survival <- frs %>%
  pull_demographic("SC\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024)
ntot <- frs %>%
  pull_demographic("NC\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024) %>%
  group_by(yr, bayes) %>%
  summarise(lci = mean(lci), mci = mean(mci), uci = mean(uci), sd = mean(sd))
lambda <- frs %>%
  pull_demographic("LAMBDA\\[") %>%
  mutate(yr = 1988 + 6*(session - 1) + yr_index) %>%
  select(yr, bayes, lci, mci, uci, sd) %>%
  filter(yr != 2024)
recruitment_betas <- frs %>%
  select(session, bayes, R_B0, R_CG, R_WM, R_WT) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "beta", values_to = "val") %>%
  group_by(session, bayes, beta) %>%
  summarize(
    lci = quantile(val, 0.025, na.rm = T),
    mci = quantile(val, 0.5, na.rm = T),
    uci = quantile(val, 0.975, na.rm = T),
    sd = sd(val, na.rm = T)
  ) %>%
  ungroup() 
survival_betas <- frs %>% 
  select(session, bayes, SC_B0, SC_CG, SC_WM, SC_WT) %>%
  pivot_longer(cols = 3:ncol(.), names_to = "beta", values_to = "val") %>%
  group_by(session, bayes, beta) %>%
  summarize(
    lci = quantile(val, 0.025, na.rm = T),
    mci = quantile(val, 0.5, na.rm = T),
    uci = quantile(val, 0.975, na.rm = T),
    sd = sd(val, na.rm = T)
  ) %>%
  ungroup() 

lambda %>% 
  select(yr, bayes, lci, mci, uci, sd) %>%
  pivot_wider(names_from = bayes, values_from = c(lci, mci, uci, sd)) %>%
  ggplot(aes(x = mci_pb, y = mci_fb, color = sd_pb, size = sd_fb)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

survival_betas %>%
  filter(session < 5) %>%
  mutate(sb = paste(session, beta)) %>%
  ggplot(aes(x = sb, y = mci, color = bayes, shape = beta)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge2(0.5))

ca_survival %>%
  filter(yr < 2013) %>%
  pivot_wider(names_from = bayes, values_from = c(lci, mci, uci, sd)) %>%
  mutate(sd_diff = sd_pb - sd_fb) %>%
  ggplot(aes(x = mci_pb, y = mci_fb, color = sd_diff)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

nt <- ntot %>%
  ggplot(aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge2(0.5)) +
  scale_color_discrete(
    name = NULL,
    labels = c("Fully integrated", "Partial Bayes")
  ) +
  xlab("Year") + ylab("Abundance") +
  theme_classic() +
  theme(legend.position = "bottom")
fs <- af_survival %>%
  ggplot(aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge2(0.5)) +
  scale_color_discrete(
    name = NULL,
    labels = c("Fully integrated", "Partial Bayes")
  ) +
  xlab(NULL) + ylab("Female Survival") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
cs <- ca_survival %>%
  ggplot(aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge2(0.5)) +
  scale_color_discrete(
    name = NULL,
    labels = c("Fully integrated", "Partial Bayes")
  ) +
  xlab(NULL) + ylab("Calf Survival") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
re <- recruitment %>%
  ggplot(aes(x = yr, y = mci, color = bayes)) +
  geom_pointrange(aes(ymin = lci, ymax = uci), position = position_dodge2(0.5)) +
  scale_color_discrete(
    name = NULL,
    labels = c("Fully integrated", "Partial Bayes")
  ) +
  xlab(NULL) + ylab("Recruitment") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
plot_grid(re, cs, fs, nt, nrow = 4, ncol = 1, rel_heights = c(3,3,3,4))
ggsave("figures//fb-pb_dem_comp.png", dpi = 600, units = "cm", height = 28, 
       width = 18)

sb_wide <- survival_betas %>%
  filter(session < 5) %>%
  select(session, beta, mci, sd, bayes) %>%
  pivot_wider(names_from = bayes, values_from = c(mci, sd)) %>%
  filter(beta != "SC_B0") %>%
  mutate(Parameter = "Calf survival") %>%
  mutate(Covariate = case_when(
    beta == "SC_CG" ~ "Puma density",
    beta == "SC_WM" ~ "SPEI (t-1)",
    beta == "SC_WT" ~ "SPEI (t)"
  ))
rb_wide <- recruitment_betas %>%
  select(session, beta, mci, sd, bayes) %>%
  pivot_wider(names_from = bayes, values_from = c(mci, sd)) %>%
  filter(beta != "R_B0") %>%
  mutate(Parameter = "Recruitment") %>%
  mutate(Covariate = case_when(
    beta == "R_CG" ~ "Puma density",
    beta == "R_WM" ~ "SPEI (t-1)",
    beta == "R_WT" ~ "SPEI (t)"
  ))

bvmed <- bind_rows(sb_wide, rb_wide) %>%
  ggplot(aes(x = mci_pb, y = mci_fb, color = Covariate, shape = Parameter)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "#555555", alpha = 0.5) +
  theme_classic() +
  xlab("Partial Bayes") + ylab("Fully Integrated") +
  labs(title = "Beta value median")
bvsd <- bind_rows(sb_wide, rb_wide) %>%
  ggplot(aes(x = sd_pb^2, y = sd_fb^2, color = Covariate, shape = Parameter)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "#555555", alpha = 0.5) +
  theme_classic() +
  xlab("Partial Bayes") + ylab("Fully Integrated") +
  labs(title = "Beta value variance") +
  theme(legend.position = "none")
plot_grid(bvmed, bvsd, nrow = 1, ncol = 2, rel_widths = c(4,3))
ggsave("figures//pb-fb_cov_comparison.png", dpi = 600, units = "cm", 
       width = 18, height = 10)
