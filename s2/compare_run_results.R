require(tidyverse); require(rjags); require(cowplot)
pull_cl_cov_name <- function(x){
  y <- x[3] %>%
    strsplit(".", fixed = T)
  return(y[[1]][1])
}
mcmc_to_tibble <- function(x){
  x %>%
    map(as_tibble) %>%
    bind_rows() %>%
    return()
}
add_cl_cov_column <- function(x, cl_cov_name){
  x %>%
    mutate(cl_cov = cl_cov_name) %>%
    return()
}
rs_file_names <- paste0("s2//results//", list.files("s2//results"))[c(4,12)]
cl_covariate_names <- strsplit(rs_file_names, "_") %>%
  map(pull_cl_cov_name) %>%
  unlist()
pu_covariate_names <- c("Logistic Growth", "Estimate Mean")
all_frs <- map(rs_file_names, readRDS) %>%
  set_names(cl_covariate_names) %>%
  map(mcmc_to_tibble) %>%
  map2(.x = ., .y = pu_covariate_names, add_cl_cov_column) %>%
  bind_rows()

covariate_summaries <- all_frs %>%
  select(cl_cov, S_cg, S_dd, S_wm, S_wt, R_cg, R_dd, R_wm, R_wt) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Covariate", values_to = "val") %>%
  group_by(cl_cov, Covariate) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95),
    p_gt0 = mean(val > 0)
  ) %>%
  ungroup() %>%
  mutate(parameter = case_when(
    Covariate %in% c("R_cg", "R_dd", "R_wm", "R_wt") ~ "recruitment",
    T ~ "survival"
  ))

# Plot Recruitment Covariates
recruitment_covariates <- covariate_summaries %>%
  filter(parameter == "recruitment") %>%
  filter(cl_cov != "null") %>%
  ggplot(aes(x = cl_cov, y = median, color = Covariate)) +
  geom_pointrange(
    aes(ymax = ucri, ymin = lcri),
    position = position_dodge2(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "33") +
  theme_classic() +
  xlab("Model") + ylab("Estimate") + labs(title = "Recruitment covariates") +
  scale_color_discrete(
    labels = c("Puma density", "Elk density", "SPEI (t-1)", "SPEI (t)")
  )

survival_covariates <- covariate_summaries %>%
  filter(parameter == "survival") %>%
  filter(cl_cov != "null") %>%
  ggplot(aes(x = cl_cov, y = median, color = Covariate)) +
  geom_pointrange(
    aes(ymax = ucri, ymin = lcri),
    position = position_dodge2(width = 0.3)) +
  geom_hline(yintercept = 0, linetype = "33") +
  theme_classic() +
  xlab("Model") + ylab("Estimate") + labs(title = "Calf survival covariates") +
  scale_color_discrete(
    labels = c("Puma density", "Elk density", "SPEI (t-1)", "SPEI (t)")
  )
plot_grid(recruitment_covariates, survival_covariates, nrow = 2, ncol = 1)
ggsave("figures//puma_covariate_comparison.png", dpi = 600, units = "cm", width = 18,
       height = 12)
recruitment_summaries <- all_frs %>%
  select(cl_cov, grep("R\\[", names(.))) %>%
  setNames(c("cl_cov", 1989:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()
survival_ca_summaries <- all_frs %>%
  select(cl_cov, grep("survival_ca\\[", names(.))) %>%
  setNames(c("cl_cov", 1989:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()
survival_af_summaries <- all_frs %>%
  select(cl_cov, grep("survival_af", names(.))) %>%
  setNames(c("cl_cov", 1989:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()
survival_am_summaries <- all_frs %>%
  select(cl_cov, grep("survival_am", names(.))) %>%
  setNames(c("cl_cov", 1989:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()
lambda_summaries <- all_frs %>%
  select(cl_cov, grep("lambda", names(.))) %>%
  setNames(c("cl_cov", 1989:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()
N_summaries <- all_frs %>%
  select(cl_cov, grep("N_tot", names(.))) %>%
  setNames(c("cl_cov", 1988:2023)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "yr", values_to = "val") %>%
  mutate(yr = as.numeric(yr)) %>%
  group_by(cl_cov, yr) %>%
  summarize(
    lcri = quantile(val, 0.05),
    median = quantile(val, 0.5),
    ucri = quantile(val, 0.95)
  ) %>%
  ungroup()

cs <- ggplot(survival_ca_summaries, aes(x = yr, y = median, color = cl_cov)) +
  geom_line() +
  xlab(NULL) + ylab("Calf Survival") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
fs <- ggplot(survival_af_summaries, aes(x = yr, y = median, color = cl_cov)) +
  geom_line() +
  xlab(NULL) + ylab("Female Survival") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
re <- ggplot(recruitment_summaries, aes(x = yr, y = median, color = cl_cov)) +
  geom_line() +
  xlab(NULL) + ylab("Recruitment") +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank())
nt <- ggplot(N_summaries, aes(x = yr, y = median, color = cl_cov)) +
  geom_line() +
  xlab("Year") + ylab("Abundance") +
  theme_classic() +
  scale_color_discrete(
    labels = c("NDVI", "NULL", "PDSI", "Precipitation", "SPEI 12m", "SPEI 3m",
               "SPEI 6m", "Temperature"),
    name = "Covariate"
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2))
plot_grid(cs, fs, re, nt, nrow = 4, ncol = 1, rel_heights = c(2,2,2,3))
ggsave("figures//covariate_demographic_comparison.png", width = 18, height = 18,
       units = "cm", dpi = 600)
