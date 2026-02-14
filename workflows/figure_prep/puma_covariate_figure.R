require(tidyverse); require(rjags)
cvdt <- readRDS("data//unscaled_covariates.rds") %>%
  select(yr, pd_reconstruction, pd_odfw_est, pd_logistic, pd_full_mean, pd_mortalities)
mortdt <- read_csv("data//deprecated//cougars//starkey_unit_mort_data.csv")
norm_rslt <- readRDS("results//fbipm_rslt_06feb2026_best_norm.rds") %>%
  map(as.matrix) %>% map(as_tibble) %>% bind_rows() %>% 
  select(grep("cdrng", names(.)))
names(norm_rslt) <- as.character(1988:2023)
get_mnsd <- function(dt, param){
  x <- dt %>%
    filter(yr > 1993 & yr < 2013) %>%
    pull(param)
  out <- c(mn = mean(x), sd = sd(x))
  return(out)
}
rescale_column <- function(dt, param){
  mnsd <- get_mnsd(dt, param)
  x <- dt %>% pull(param)
  out <- (x - mnsd[1]) / mnsd[2]
  return(out)
}

puma <- cvdt %>%
  mutate(
    rcns = rescale_column(., "pd_reconstruction"),
    lgst = rescale_column(., "pd_logistic"),
    odfw = rescale_column(., "pd_odfw_est"),
    mort = rescale_column(., "pd_mortalities")
  ) %>%
  select(yr, rcns, lgst, odfw, mort) %>%
  mutate(mnes = rowMeans(select(., 2:4), na.rm = T)) %>%
  filter(yr > 1987) %>%
  pivot_longer(2:ncol(.), names_to = "prm", values_to = "val") %>%
  mutate(parameter = case_when(
    prm == "rcns" ~ "Reconstruction",
    prm == "lgst" ~ "Logistic Growth",
    prm == "odfw" ~ "ODFW Estimate",
    prm == "mnes" ~ "Mean Index",
    prm == "mort" ~ "Mortalities"
  ))
puma %>%
  filter(prm %in% c("rcns", "odfw", "mort")) %>%
  ggplot(aes(x = yr, y = val, color = parameter, linetype = parameter)) +
  # geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.title = element_blank()) +
  xlab("Year") + ylab("Scaled Value") +
  theme(legend.position = "bottom")
ggsave(
  "figures//puma_data_sources.png",
  width = 14,
  height = 7,
  units = "cm",
  dpi = 600
)
mortdt %>% 
  filter(Age < 20) %>%
  ggplot(aes(Age)) +
  geom_histogram(binwidth = 1, color = "#000", fill = "#999") +
  theme_classic() + ylab("Count") + xlab("Age at Death")
ggsave(
  "figures//puma_mortality_ages.png",
  width = 8,
  height = 8,
  units = "cm",
  dpi = 600
)

lngdt <- norm_rslt %>%
  mutate(stpn = 1:nrow(.)) %>%
  pivot_longer(1:(ncol(.)-1), names_to = "yr", values_to = "val")
rngmn <- lngdt %>% filter(yr > 1993 & yr < 2013) %>% pull(val) %>% mean()
rngsd <- lngdt %>% filter(yr > 1993 & yr < 2013) %>% pull(val) %>% sd()

cgnorm <- lngdt %>%
  mutate(y = (val - rngmn) / rngsd) %>%
  group_by(yr) %>%
  summarise(
    lci = quantile(y, 0.025),
    val = quantile(y, 0.5),
    uci = quantile(y, 0.975)
  ) %>%
  mutate(yr = as.numeric(yr))
cgnormadd <- cgnorm %>%
  mutate(prm = "norm", parameter = "Random Normal") %>%
  select(yr, prm, val, parameter)
puma_plot <- puma %>%
  bind_rows(cgnormadd) %>%
  filter(prm %in% c("rcns", "odfw", "lgst", "mnes", "norm")) %>%
  ggplot(aes(x = yr, y = val, color = parameter, linetype = parameter)) +
  # geom_point() +
  geom_line() +
  geom_ribbon(
    data = cgnorm, 
    aes(ymin = lci, ymax = uci), 
    alpha = 0.2,
    color = "#ccc",
    fill = "#E69F00",
    linetype = 1,
    linewidth = 0.1
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "Reconstruction" = "#800080",
    "Random Normal" = "#E69F00",
    "ODFW Estimate" = "#009E73",
    "Mean Index" = "#D55E00",
    "Logistic Growth" = "#0072B2"
  )) +
  theme(legend.title = element_blank()) +
  xlab("Year") + ylab("Scaled Value") +
  theme(legend.position = "bottom") 
ggsave(
  "figures//puma_covariate_values.png",
  width = 18,
  height = 9,
  units = "cm",
  dpi = 600
)
